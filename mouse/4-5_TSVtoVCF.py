import argparse
import textwrap
import re

# Parse arguments
parser = argparse.ArgumentParser(description=textwrap.dedent("""\
                                    Convert files from SComatic TSV format to VCFv4.3

                                    ---------------------------------------------------------
                                    IMPORTANT
                                    - This script has been tested only with output from
                                    SComatic downloaded from GitHub commit f515f4e
                                    ---------------------------------------------------------
                                    """), formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--add-info', action='store_true', help="Add the extra columns from SComatic TSV as key-value pairs in the VCF INFO column")
parser.add_argument('-c', '--add-celltypes', action='store_true', help="Add a FORMAT column and cell type information as column samples to the VCF file")
parser.add_argument('file', action='store', help="SComatic TSV pre-filtered to include only PASS variants")
parser.add_argument('output', action='store', help="VCF file to write")

args = parser.parse_args()

# MAIN
vcf_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

with open(args.file, 'r') as infile, open(args.output, 'w') as outfile:
    
    # Manually add the meta-information in the header
    meta_info = textwrap.dedent("""\
            ##fileformat=VCFv4.3
            ##convertedFrom=SComaticTSV
            ##conversionTool=scomatic-to-vcf.py
            ##reference=~/bin/SComatic-main/reference_genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa
            ##contig=<ID=1,length=195471971,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=2,length=182113224,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=3,length=160039680,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=4,length=156508116,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=5,length=151834684,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=6,length=149736546,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=7,length=145441459,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=8,length=129401213,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=9,length=124595110,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=10,length=130694993,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=11,length=122082543,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=12,length=120129022,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=13,length=120421639,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=14,length=124902244,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=15,length=104043685,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=16,length=98207768,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=17,length=94987271,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=18,length=90702639,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=19,length=61431566,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=22,length=61431566,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=X,length=171031299,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=Y,length=91744698,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=MT,length=16299,assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=GL456370.1,length="",assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##contig=<ID=JH584304.1,length="",assembly=GRCm38,species="Mus musculus",taxonomy=x>
            ##FILTER=<ID=Multiple_cell_types,Description="Quality below 10">
            ##FILTER=<ID=Cell_type_noise,Description="Quality below 10">
            ##FILTER=<ID=Min_cell_types,Description="Quality below 10">
            ##FILTER=<ID=Noisy_site,Description="Quality below 10">
            ##FILTER=<ID=Clustered,Description="Quality below 10">
            ##FILTER=<ID=LC_Upstream,Description="Quality below 10">
            ##FILTER=<ID=Multi-allelic,Description="Quality below 10">
            ##FILTER=<ID=LC_Downstream,Description="Quality below 10">
            ##INFO=<ID=Cell_types,Number=.,Type=String,Description="Cell type/s with the variant">
            ##INFO=<ID=Up_context,Number=.,Type=String,Description="Up-stream bases in reference (4 bases)">
            ##INFO=<ID=Down_context,Number=.,Type=String,Description="Down-stream bases in reference (4 bases)">
            ##INFO=<ID=N_ALT,Number=1,Type=Integer,Description="Cell type/s with the variant">
            ##INFO=<ID=Dp,Number=.,Type=Integer,Description="Depth of coverage (reads) in the cell type supporting the variant">
            ##INFO=<ID=Nc,Number=.,Type=Integer,Description="Number of distinct cells found in the cell type with the mutation">
            ##INFO=<ID=Bc,Number=.,Type=Integer,Description="Number of reads (base count) supporting the variants in the cell type with the mutation">
            ##INFO=<ID=Cc,Number=.,Type=Integer,Description="Number of distinct cells supporting the variant in the cell type with the mutation">
            ##INFO=<ID=VAF,Number=.,Type=Float,Description="Variant allele frequency of variant in the cell type with the mutation">
            ##INFO=<ID=CCF,Number=.,Type=Float,Description="Cancer cell fraction (fraction of ditinct cells) supporting the alternative allele in the cell type with the mutation">
            ##INFO=<ID=BCp,Number=.,Type=Float,Description="Beta-binomial p-value for the variant allele (considering read counts)">
            ##INFO=<ID=CCp,Number=.,Type=Float,Description="Beta-binomial p-value for the variant allele (considering cell counts)">
            ##INFO=<ID=Cell_types_min_BC,Number=1,Type=Integer,Description="Number of cell types with a minimum number of reads covering a site">
            ##INFO=<ID=Cell_types_min_CC,Number=1,Type=Integer,Description="Number of cell types with a minimum number of distinct cells found in a specific site">
            ##INFO=<ID=Rest_BC,Number=3,Type=Float,Description="Base counts (reads) supporting other alternative alleles in this site. BC;DP;P-value (betabin)">
            ##INFO=<ID=Rest_CC,Number=3,Type=Float,Description="Cell counts supporting other alternative alleles in this site. CC;NC;P-value (betabin)">
            ##INFO=<ID=Fisher_p,Number=.,Type=Float,Description="Strand bias test. Fisher exact test p-value between forward and reverse reads in variant and reference allele">
            ##INFO=<ID=Cell_type_Filter,Number=.,Type=String,Description="Filter status of the variant site in each cell type">
            ##INFO=<ID=Multiple_cell_types,Number=.,Type=String,Description="Filter status of the variant site in each cell type">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of coverage">
            ##FORMAT=<ID=NC,Number=1,Type=Integer,Description="Number of different cells">
            ##FORMAT=<ID=CC,Number=1,Type=String,Description="Cell counts [A|C|T|G|I|D|N|O], where D means deletion, I insertion and O other type of character">
            ##FORMAT=<ID=BC,Number=1,Type=String,Description="Base counts [A|C|T|G|I|D|N|O], where D means deletion, I insertion and O other type of character">
            ##FORMAT=<ID=BCf,Number=1,Type=String,Description="Base counts in forward reads [A|C|T|G|I|D|N|O], where D means deletion, I insertion and O other type of character">
            ##FORMAT=<ID=BCr,Number=1,Type=String,Description="Base counts in reverse reads [A|C|T|G|I|D|N|O], where D means deletion, I insertion and O other type of character">
            ##FORMAT=<ID=BQ,Number=1,Type=String,Description="No info">
            """)

    outfile.writelines(meta_info)

    # Add header and data lines with variant information
    for line in infile:
        if not line.startswith("##"):
            if line.startswith("#"):
                # Parse header
                header = line.strip('#\n').split('\t')

                # Write VCF header
                if args.add_celltypes:
                    # Add header with "FORMAT" column and the rest of specific cell type columns
                    # Cell type columns appear after column 24 (using 0-based indexing)
                    vcf_cols.extend(header[25:])
                    outfile.write("#{}\n".format('\t'.join(vcf_cols)))
                else:
                    # Without "FORMAT" column
                    outfile.write("#{}\n".format('\t'.join(vcf_cols[:8])))
            else:
                # Save each SComatic record to its corresponding VCF column
                fields = line.rstrip().split('\t')

                col_chrom = fields[0]
                col_pos = fields[1]
                col_id = "."
                col_ref = fields[3]
                col_alt = ",".join(set(fields[4].split(',')))  # VCF requires no duplicated ALT values
                col_qual = "."
                col_filter = re.sub(",", ";", fields[5])

                # Populate info column
                if args.add_info:
                    # Add all columns with their associated heading between column 6 (FILTER) and 24 (INFO) (0-based index) 
                    info_list = [ "{}={}".format(col, re.sub(";", ",", field)) for col, field in zip(header[6:24], fields[6:24]) ]
                    col_info = ";".join(info_list)
                else:
                    col_info = "."
                
                # Populate sample columns
                if args.add_celltypes:
                    # Get INFO column from SComatic output. This will be FORMAT column in VCF.
                    col_format = re.sub("[|]", ":", fields[24])
                    # Add each cell type column as a sample. From column 25 to the last column
                    # Make some changes to comply with VCF 4.3 scpecification to cell type information
                    col_samples = '\t'.join(fields[25:])
                    t = col_samples.maketrans({'|': ':', ':': '|'})  # interchange ":" and "|" because of the VCF specification
                    col_samples = col_samples.translate(t)
                    col_samples = re.sub("NA", re.sub("[A-Za-z]{2,3}", ".", col_format), col_samples)  # chage the NA values for the MISSING value "." of VCF 

                newfields = [col_chrom, col_pos, col_id, col_ref, col_alt, col_qual, col_filter, col_info, col_format, col_samples]

                newline = "\t".join(newfields)
                
                outfile.write(newline + "\n")
