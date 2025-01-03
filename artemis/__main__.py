import gzip
import re
import sys
import pysam
import argparse
import os

from pybedtools import BedTool
import logging

parser = argparse.ArgumentParser(prog="artemis", usage="artemis -h", description="Run the artemis pipeline", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
handler = logging.StreamHandler()
handler.terminator = ""

logging.basicConfig(level=logging.INFO, format='%(message)s', handlers=[handler])

hg={}
def loadhg(args):

    logging.info("\nLoading genome %s... "%args.reference.name)

    f=args.reference.read()

    if f[:2] == b'\x1f\x8b': #gzip magic number
        logging.info('Reading gzipped genome file...')
        genome=gzip.decompress(f).decode('utf-8')
    else:
        genome=f.decode('utf-8')

    nseq=0

    seq=""
    for line in genome.split("\n"):
        if len(line)>0:
            if line[0]=='>':
                if seq!="":
                    hg[ctgname]=seq
                    nseq+=len(seq)
                ctgname=line[1:].split()[0]

                seq=""
            else:
                seq+=line[:-1]

    if seq!="":
        hg[ctgname]=seq
        nseq+=len(seq)

    logging.info(f" done ({nseq}bp).\n")

revcomptable = str.maketrans("acgtACGTRY","tgcaTGCAYR")
def revcomp(s):
    return s.translate(revcomptable)[::-1]

def intersect_seeds_with_vcf(args):    
    pamsites = []
    seedsites= []

    logging.info("Scanning for PAM sites in genome... ")
    p=re.compile("(TTT[A,C,G])|([C,G,T]AAA)")
    for chrom in hg.keys():
        logging.info(f"{chrom} ")
        for m in p.finditer(hg[chrom]):
            if chrom=='MT':
                chrom='M'
            o='-' if m.groups()[0]==None else '+'
            chrom, start, end, orient = chrom, m.start(), m.end(), o

            pamsites.append((chrom, start, end, orient))
    
    logging.info(f"done ({len(pamsites)} pamsites).\n")

    pams = BedTool(pamsites)
    pams.saveas("allpamsites.bed")

    if args.kg!=None:
        logging.info("Intersecting PAM sites with 1000 genomes data... ")
        kg=BedTool(open(args.kg.name))
        pams = pams.intersect(kg,  wa=True, wb=False, v=True)
        pams.saveas("filtpamsites.bed")
        logging.info("done.\n")

    logging.info("Sorting and merging %d PAM associated seed sites... "%len(pamsites))
    for p in pams:
        chrom, start, end, orient = p
        start, end = int(start), int(end)
        if orient == "+":
            seedsites.append((chrom, end, end + args.seedsize, orient))
        else:
            seedsites.append((chrom, start-args.seedsize, start, orient))
    seeds = BedTool(seedsites)
    seeds = seeds.sort().merge()
    logging.info("done.\n")

    logging.info("Storing seed regions... ")
    seeds.saveas("seedsites.bed")
    logging.info("done.\n")

    clinvarvcf=open(args.vcf.name).read()

    logging.info("Intersecting seed regions with input VCF file... ")
    inputvcf=BedTool(clinvarvcf, from_string=True)
    intersected_bed = inputvcf.intersect(seeds, wa=True, wb=False, u=True)
    logging.info("done.\n")

    logging.info("Writing all CAS12 seed intersecting variants to new VCF file with annotations... ")

    header=""
    num_input_snvs=0
    for l in clinvarvcf.split("\n"):
        if len(l)>0:
            if l[0]!='#':
                num_input_snvs+=1
            else:
                header+=l+"\n"

    logging.info("Number of SNVs in input VCF: %d\n"%num_input_snvs)

    intersected_bed.saveas(".intersection.vcf", trackline=header)

    if args.outputfile == None:
        args.outputfile = sys.stdout

    # Write all CAS12 intersecting variants to a new annotated VCF file
    with pysam.VariantFile(".intersection.vcf") as vcfreader:

        posPAM =  ["TTTC", "TTTA", "TTTG"]
        negPAM  = ["CAAA", "TAAA", "GAAA"]
        pamsize=4
        w=args.seedsize+pamsize

        n=0

        vcfreader.header.info.add('SNP_position', 1, 'String', 'Position of SNP relative to PAM site')
        vcfreader.header.info.add('CRISPR_ref_target_sequence', 1, 'String', 'CRISPR target sequence to address the reference (WT) allele.')
        vcfreader.header.info.add('CRISPR_alt_target_sequence', 1, 'String', 'CRISPR target sequence to address the alternative (pathogenic) allele.')
        vcfreader.header.info.add('STRAND', 1, 'String', 'Whether PAM site is found on the positive (+) or negative (-) strand.')
        
        with pysam.VariantFile(args.outputfile,'w',header=vcfreader.header) as vcfwriter:
            
            for rec in vcfreader:
                
                if rec.ref!=None and rec.alts!=None and rec.pos!=None and rec.chrom!=None:

                    if len(rec.ref) == 1 and len(rec.alts)==1 and len(rec.alts[0]) == 1: #only SNVs!

                        chrom=rec.chrom
                        chrom_nochr=chrom.replace('chr','').replace('M','MT')
                        pos=rec.pos
                        ref=rec.ref
                        alt=rec.alts[0]
                        
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
                        n+=1
    logging.info("done.\n")

    logging.info("%.3f%% variants in the inputvcf (%d/%d) intersect a CAS12 seed site.\n"%(float(n/num_input_snvs)*100,n,num_input_snvs))
    os.remove(".intersection.vcf")

def main():
    parser.add_argument(dest="reference", type=argparse.FileType('rb'), help="Path to the human reference genome file")
    parser.add_argument(dest="vcf", type=argparse.FileType('r'), help="Path to the a vcf file, e.g. a clinvar vcf file")

    parser.add_argument("--exclude", dest="kg", type=argparse.FileType('r'), help="Exclude PAMs that intersect with variants in this vcf (e.g. 1000 genomes data)", required=False, default=None)
    parser.add_argument("-o", dest="outputfile", type=argparse.FileType('w'), help="Where to write annotated output vcf file (default: stdout)", required=False, default=None)
    parser.add_argument("--seedsize", type=int, default=5, help="Specify the size of the seed region")

    args = parser.parse_args()

    loadhg(args)

    intersect_seeds_with_vcf(args)

if __name__ == '__main__':
    main()
