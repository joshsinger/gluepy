../bin/neo4j-admin import \
     --nodes=Sequence.csv \
     --nodes=Snp.csv \
     --nodes=AaReplacement.csv \
     --relationships=CONTAINS_SNP.csv \
     --relationships=CONTAINS_AA_REP.csv \
     --trim-strings=true