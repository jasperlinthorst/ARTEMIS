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
                seq+=line#.upper() #consider uppercasing the genome to include pams in softmasked regions...

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
    p=re.compile(args.pamregex)
    for chrom in hg.keys():
        logging.info(f"{chrom} ")
        for m in p.finditer(hg[chrom]):

            o='-' if m.groups()[0]==None else '+'
            start=m.start()
            end=m.start() + len(m.group(1)) if m.group(1)!=None else m.start()+len(m.group(2))
            
            # print(chrom,start,end,hg[chrom][start:end])
            # print(m.groups(), start, end)
            
            # if o=='+':
            #     assert(hg[chrom][start:end] in ["TTTC", "TTTA", "TTTG"])
            # else:
            #     assert(hg[chrom][start:end] in ["CAAA", "TAAA", "GAAA"])
            
            pamsites.append((chrom, start, end, o))
    
    logging.info(f"done ({len(pamsites)} pamsites).\n")

    pams = BedTool(pamsites)
    pams.saveas("allpamsites.bed")

    if args.kg!=None:
        logging.info("Intersecting PAM sites with exclusion vcf... ")
        kg=BedTool(open(args.kg.name).read(),from_string=True)
        pams = pams.intersect(kg,  wa=True, wb=False, v=True)
        pams.saveas("filtpamsites.bed")
        logging.info("done.\n")
    
    logging.info("Sorting and merging %d PAM associated seed sites... "%len(pams))
    for p in pams:
        chrom, start, end, orient = p
        start, end = int(start), int(end)
        if start > args.seedsize and end < len(hg[chrom]) - args.seedsize:
            if orient == "+":
                seedsites.append((chrom, end, end + args.seedsize, orient))
            else:
                seedsites.append((chrom, start-args.seedsize, start, orient))
        else:
            logging.info(f"Skipping PAM site at {chrom}:{start}-{end} because it is too close to the edge of the chromosome.\n")
    
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

    intersected_bed.saveas(".intersection.vcf", trackline=header)

    if args.outputfile == None:
        args.outputfile = sys.stdout

    # Write all CAS12 intersecting variants to a new annotated VCF file
    with pysam.VariantFile(".intersection.vcf") as vcfreader:

        gRNAsize=25
        w=gRNAsize #should be max seedsize + PAMsize, but theoretically regex matches could be variably sized, so now we scan more than we need, but we don't need an additional parameter
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
                        pos=rec.pos
                        ref=rec.ref
                        alt=rec.alts[0]
                                                
                        region=hg[chrom][int(pos)-w:int(pos)+w-1] #determine w or use valid default
                        p=re.compile(args.pamregex)
                        
                        rpclosest=None
                        pamclosest=None
                        for m in p.finditer(region):

                            if m.groups()[0]!=None: #PAM on positive strand
                                pam=m.group(1)
                                strandorientation='+'
                                pamsize=len(pam)
                                if m.start()+len(pam)<w: #valid because its upstream of the variant
                                    rp=w-m.start() #position of the variant relative to the start of the PAM site
                                else: continue #pams that intersect the SNP are not valid (?)
                            else: # PAM on negative strand
                                pam=m.group(2)
                                strandorientation='-'
                                if m.start()>=w: #valid because its downstream of the variant
                                    rp=m.start()+1+len(pam)-w #position of the variant relative to the start of the PAM site
                                else: continue #pams that intersect the SNP are not valid (?)
                            if rpclosest==None or rp<rpclosest:
                                rpclosest=rp
                                pamclosest=pam
                                orientation=strandorientation
                                if orientation=='+':
                                    pamstart=int(pos)-rpclosest
                                    reftarget=hg[chrom][pamstart:pamstart+gRNAsize]
                                    alttarget=reftarget[:(rpclosest-1)]+alt+reftarget[rpclosest:]
                                else:
                                    pamstart=int(pos)-1+rpclosest
                                    reftarget=revcomp(hg[chrom][pamstart-gRNAsize:pamstart])
                                    alttarget=reftarget[:(rpclosest-1)]+revcomp(alt)+reftarget[rpclosest:]

                                # print(n,chrom,ref,alt,pos,"is targetable with CAS12 through PAM site: %s at position +%d!"%(pamclosest,rpclosest-len(pamclosest)),downstream,r,upstream,region)
                        
                        rec.info['SNP_position']='PAM+%d'%(rpclosest-len(pamclosest))
                        rec.info['CRISPR_ref_target_sequence']=reftarget.upper()
                        rec.info['CRISPR_alt_target_sequence']=alttarget.upper()
                        rec.info['STRAND']=orientation
                        
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
    parser.add_argument("--pamregex", type=str, default='(?=(TTT[ACG]))|(?=([CGT]AAA))', help="Regular expression to search for PAM sites in the genome, default for CAS12a: first group should match forward strand, second group the reverse complement")
    
    args = parser.parse_args()

    loadhg(args)

    intersect_seeds_with_vcf(args)

if __name__ == '__main__':
    main()
