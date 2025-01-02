import gzip
import re
import sys
import pysam
import argparse
from pybedtools import BedTool
import logging

parser = argparse.ArgumentParser(prog="artemis", usage="artemis -h", description="Run the artemis pipeline", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stderr)

hg={}
def loadhg(args):
    logging.info("\nLoading human reference genome")

    if args.reference.name.endswith(".gz"):
        fopen=gzip.open
    else:
        fopen=open
    
    with fopen(args.reference.name, "rt") as genome_file:
        seq,chrom="",None
        for line in genome_file:
            if line[0]=='>':
                if seq!="" and chrom!=None:
                    hg[chrom]=seq
                p=re.compile("(.*) Homo sapiens chromosome (\d{1,2}|X|Y), GRCh38.p14 Primary Assembly")
                match=p.search(line)
                if match!=None: #actual primary assembly
                    ctgname=match.group(1)
                    chrom=match.group(2)
                else:
                    if line.find('mitochondrion')!=-1:
                        chrom='M'
                    else:
                        chrom=None  
                seq=""
            else:
                seq+=line[:-1]
    if seq!="" and chrom!=None:
        hg[chrom]=seq

revcomptable = str.maketrans("acgtACGTRY","tgcaTGCAYR")
def revcomp(s):
    return s.translate(revcomptable)[::-1]

def intersect_seeds_with_clinvar(args):
    

    # Step 3: Create a set of positions (chrom, start, end) from the VCF file
    vcf_positions = set()

    if args.kg!=None:
        logging.info("Load 1000 genomes data")
        kgvcf = pysam.VariantFile(args.kg.name)
        for record in kgvcf.fetch():
            # Collect positions where AF is between 0.01 and 0.99
            if "AF" in record.info:
                af_values = record.info["AF"]
                if any(0.01 <= af <= 0.99 for af in af_values): #TODO: parameterize this
                    vcf_positions.add((record.chrom, record.pos - 1, record.pos))  # VCF is 1-based, BED is 0-based

    # Step 4: Perform the intersection to filter out BED intervals that overlap with VCF positions
    no_overlap_bed = []
    # for interval in bed:

    logging.info("Scan PAM sites in human genome")
    p=re.compile("(TTT[A,C,G])|([C,G,T]AAA)")
    for chrom in hg.keys():
        for m in p.finditer(hg[chrom]):
            if chrom=='MT':
                chrom='M'
            o='-' if m.groups()[0]==None else '+'
            chrom, start, end, orient = chrom, m.start(), m.end(), o

            # chrom, start, end, orient = interval.chrom, interval.start, interval.end, interval.name
            # Check if the interval overlaps with any of the VCF positions
            if not any(chrom == vcf_chrom and start < vcf_end and end > vcf_start
                    for vcf_chrom, vcf_start, vcf_end in vcf_positions):
                # if no overlap, determine the seed region
                if orient == "+":
                    no_overlap_bed.append((chrom, end, end + args.seedsize, orient,'test'))
                else:
                    no_overlap_bed.append((chrom, start-args.seedsize, start, orient,'test'))
    
    no_overlap_bed_tool = BedTool(no_overlap_bed)

    logging.info("Sort and merge %d PAM associated seed sites"%len(no_overlap_bed))
    merged_bed_tool = no_overlap_bed_tool.sort().merge()

    merged_bed_tool.saveas("seedsites.bed")

    

    # clinvarstr="".join(gzip.open(args.clinvar.name,"rt").readlines())

    clinvarvcf = BedTool(args.clinvar.name)
    
    #clinvar_bed = BedTool([(rec.chrom, rec.pos - 1, rec.pos) for rec in pysam.VariantFile(args.clinvar).fetch() if "CLNSIG" in rec.info and "Pathogenic" in rec.info["CLNSIG"] and 'CLNVC' in rec.info and rec.info["CLNVC"] == 'single_nucleotide_variant'])

    intersected_bed = clinvarvcf.intersect(merged_bed_tool, wa=True, wb=False, u=True)

    intersected_bed.saveas("clinvar_intersect.bed")

    
    logging.info("Writing CRISPR gRNA annotated clinvar VCF file for intersecting clinvar pathogenic SNVs")

    with pysam.VariantFile(args.clinvar.name, "r") as clinvarvcf:
        header=str(clinvarvcf.header)

    intersected_bed.saveas("clinvar_intersection.vcf", trackline=header)

    # Write all pathogenic SNV variants in clinvar that intersect with a CAS12 seed to a new VCF file with specific annotations
    with pysam.VariantFile("clinvar_intersection.vcf") as vcfreader:

        posPAM =  ["TTTC", "TTTA", "TTTG"]
        negPAM  = ["CAAA", "TAAA", "GAAA"]
        pamsize=4
        w=args.seedsize+pamsize

        n=0

        vcfreader.header.info.add('SNP_position', 1, 'String', 'Position of SNP relative to PAM site')
        vcfreader.header.info.add('CRISPR_ref_target_sequence', 1, 'String', 'CRISPR target sequence to address the reference (WT) allele.')
        vcfreader.header.info.add('CRISPR_alt_target_sequence', 1, 'String', 'CRISPR target sequence to address the alternative (pathogenic) allele.')
        vcfreader.header.info.add('STRAND', 1, 'String', 'Whether PAM site is found on the positive (+) or negative (-) strand.')
        
        with pysam.VariantFile(sys.stdout,'w',header=vcfreader.header) as vcfwriter:
            
            for rec in vcfreader:
                
                if "CLNSIG" in rec.info and "Pathogenic" in rec.info["CLNSIG"]: #Only write pathogenic variants
                    if 'CLNVC' in rec.info and rec.info["CLNVC"]=='single_nucleotide_variant': #Only write SNVs
                        # chrom, pos = rec.chrom, rec.pos - 1  # VCF is 1-based, BED is 0-based
                        # for interval in merged_bed_tool:
                            # if chrom == interval.chrom and pos >= interval.start and pos < interval.end:

                        chrom=rec.chrom
                        chrom_nochr=chrom.replace('chr','').replace('M','MT')
                        gene=rec.info['GENEINFO']
                        pos=rec.pos
                        ref=rec.ref
                        alt=rec.alts[0]
                        #af=rec.info['AF_EXAC']
                        
                        site=hg[chrom_nochr][int(pos)-w-1:int(pos)+w]
                        r=hg[chrom_nochr][int(pos)-1]
                        downstream=hg[chrom_nochr][int(pos)-w:int(pos)-1]
                        upstream=hg[chrom_nochr][int(pos):int(pos)+w-1]
                        
                        rpclosest=100
                        for pam in posPAM:
                            p=downstream.find(pam)
                            if p!=-1:
                                rp=w-pamsize-p
                                if rp<rpclosest:
                                    strandorientation='+'
                                    rpclosest=rp
                                    pamstart=int(pos)-rpclosest-4
                                    reftarget=hg[chrom_nochr][pamstart:pamstart+25]
                                    alttarget=reftarget[:4+(rpclosest-1)]+alt+reftarget[4+rpclosest:]
                                    
                                #print(n,chrom,gene,af,ref,alt,pos,"is targetable with CAS12 through PAM site: %s at position +%d!"%(pam,w-pamsize-p),downstream,r,upstream)
                        
                        for pam in negPAM:
                            p=upstream.find(pam)
                            if p!=-1:
                                rp=p+1
                                if rp<rpclosest:
                                    strandorientation='-'
                                    rpclosest=rp
                                    pamstart=int(pos)-1+rpclosest+4
                                    reftarget=revcomp(hg[chrom_nochr][pamstart-25:pamstart])
                                    alttarget=reftarget[:4+(rpclosest-1)]+revcomp(alt)+reftarget[4+rpclosest:]
                                
                                #print(n,chrom,gene,af,ref,alt,pos,"is targetable with CAS12 through PAM site: %s at position +%d!"%(pam,p+1),downstream,r,upstream)
                        
                        rec.info['SNP_position']='PAM+%d'%rpclosest
                        rec.info['CRISPR_ref_target_sequence']=reftarget.upper()
                        rec.info['CRISPR_alt_target_sequence']=alttarget.upper()
                        rec.info['STRAND']=strandorientation
                        
                        vcfwriter.write(rec)

def main():
    parser.add_argument(dest="clinvar", type=argparse.FileType('r'), help="Path to the clinvar vcf file")
    parser.add_argument(dest="reference", type=argparse.FileType('r'), help="Path to the human reference genome file")

    parser.add_argument("--1kg", dest="kg", type=argparse.FileType('r'), help="Path to 1000 genomes vcf file", required=False, default=None)
    parser.add_argument("--seedsize", type=int, default=5, help="Specify the size of the seed region")

    args = parser.parse_args()

    loadhg(args)

    intersect_seeds_with_clinvar(args)

if __name__ == '__main__':
    main()
