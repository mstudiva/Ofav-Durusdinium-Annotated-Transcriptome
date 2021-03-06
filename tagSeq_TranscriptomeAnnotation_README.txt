# Transcriptome Annotation, version November 4, 2021
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on FAU's HPC (KoKo)
# for use in generating transcriptome annotation files for Orbicella faveolata
# also includes the separation of reads associated with O. faveolata and Durusdinium transcriptomes


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- email@gmail.com by your actual email;
#	- username with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
# setup

# To install Bioperl in your bin directory, please follow these instructions:
cd bin
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/z0on/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes
rm launcher_creator.py

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/mstudiva/Ofav-Durusdinium-Annotated-Transcriptome.git
mv Ofav-Durusdinium-Annotated-Transcriptome/* .
rm -rf Ofav-Durusdinium-Annotated-Transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
# getting transcriptomes

# O. faveolata transcriptome (April 2015)
wget https://www.dropbox.com/s/0o5ntnlymyyhzkp/Ofaveolata_transcriptome.fasta
mv Ofaveolata_transcriptome.fasta Ofaveolata.fasta

# Durusdinium transcriptome (November 2020)
wget https://marinegenomics.oist.jp/symbd/download/102_symbd_transcriptome_nucl.fa.gz
gzip -d 102_symbd_transcriptome_nucl.fa.gz
mv 102_symbd_transcriptome_nucl.fa Durusdinium.fasta

# use the stream editor to find and replace all instances of "comp" with "Ofaveolata" in the host transcriptome, and "TRINITY" with "Durusdinium" in the symbiont transcriptome
sed -i 's/comp/Ofaveolata/g' Ofaveolata.fasta
sed -i 's/TRINITY_DN/Durusdinium/g' Durusdinium.fasta

# transcriptome statistics
echo "seq_stats.pl Ofaveolata.fasta > seqstats_Ofaveolata.txt" > seq_stats
echo "seq_stats.pl Durusdinium.fasta > seqstats_Durusdinium.txt" >> seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch seq_stats.slurm

nano seqstats_Ofaveolata.txt

Ofaveolata.fasta
-------------------------
178943 sequences.
1100 average length.
38110 maximum length.
201 minimum length.
N50 = 2218
196.8 Mb altogether (196757464 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

nano seqstats_Durusdinium.txt

Durusdinium.fasta
-------------------------
122923 sequences.
1097 average length.
11620 maximum length.
201 minimum length.
N50 = 1624
134.8 Mb altogether (134792983 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
# uniprot annotations with blast

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# unzipping
gunzip uniprot_sprot.fasta.gz &

# indexing the fasta database
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch mdb.slurm

# splitting the transcriptome into 200 chunks
splitFasta.pl Ofaveolata.fasta 200
splitFasta.pl Durusdinium.fasta 200

# blasting all 200 chunks to uniprot in parallel, 4 cores per chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 6:00:00 -q shortq7 -e email@gmail.com
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l
# you should end up with the same number of queries as sequences from the seq_stats script (301866 sequences)

# combining all blast results
cat subset*br > myblast.br
mv subset* ~/annotate/backup/

# for trinity-assembled transcriptomes: annotating with "Durusdinium" or "Ofaveolata" depending on if component is from symbiont or host (=component)
grep ">" Ofaveolata.fasta | perl -pe 's/>Ofaveolata(\d+)(\S+)\s.+/Ofaveolata$1$2\tOfaveolata$1/'>Ofaveolata_seq2iso.tab
cat Ofaveolata.fasta | perl -pe 's/>Ofaveolata(\d+)(\S+).+/>Ofaveolata$1$2 gene=Ofaveolata$1/'>Ofaveolata_iso.fasta

grep ">" Durusdinium.fasta | perl -pe 's/>Durusdinium(\d+)(\S+)\s.+/Durusdinium$1$2\tDurusdinium$1/' > Durusdinium_seq2iso.tab
cat Durusdinium.fasta | perl -pe 's/>Durusdinium(\d+)(\S+).+/>Durusdinium$1$2 gene=Durusdinium$1/' > Durusdinium_iso.fasta


#-------------------------
# extracting coding sequences and corresponding protein translations:
echo "perl ~/bin/CDS_extractor_v2.pl Ofaveolata_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e email@gmail.com
sbatch cddd

echo "perl ~/bin/CDS_extractor_v2.pl Durusdinium_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e email@gmail.com
sbatch cddd


#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# selecting the longest contig per isogroup (also renames using isogroups based on Ofaveolata and Durusdinium annotations):
fasta2SBH_Ofav.pl Ofaveolata_iso_PRO.fas >Ofaveolata_out_PRO.fas

fasta2SBH_Ofav.pl Durusdinium_iso_PRO.fas >Durusdinium_out_PRO.fas

# scp your *_out_PRO.fas file to laptop, submit it to
http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp username@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*_out_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# Ofav status: go on web to http://eggnog-mapper.embl.de/job_status?jobname=MM_hsd55r5q
# symD status: go on web to http://eggnog-mapper.embl.de/job_status?jobname=MM_12stw4y1

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_hsd55r5q/out.emapper.annotations
wget http://eggnog-mapper.embl.de/MM_12stw4y1/out.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Ofaveolata_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Ofaveolata_iso2geneName.tab

awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Durusdinium_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Durusdinium_iso2geneName.tab


#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Ofaveolata_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Ofaveolata_iso2kogClass1.tab > Ofaveolata_iso2kogClass.tab

awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Durusdinium_iso2kogClass1.tab
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Durusdinium_iso2kogClass1.tab > Durusdinium_iso2kogClass.tab


#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
srun fasta2SBH_Ofav.pl Ofaveolata_iso.fasta >Ofaveolata_4kegg.fasta

srun fasta2SBH_Ofav.pl Durusdinium_iso.fasta >Durusdinium_4kegg.fasta

# scp *4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*4kegg.fasta .
# use web browser to submit Ofaveolata_Durusdinium_4kegg.fasta file to KEGG's KAAS server ( http://www.genome.jp/kegg/kaas/ )
# select SBH method, upload nucleotide query
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1636053354&key=MYtZ30dV
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1636056150&key=xK3WN2G8

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1636053354/query.ko
wget https://www.genome.jp/tools/kaas/files/dl/1636056150/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Ofaveolata_iso2kegg.tab

cat query.ko | awk '{if ($2!="") print }' > Durusdinium_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
# file transfer

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/* .
