# BIP


BIP_DATA="$(pwd)" docker compose -f compose.yaml   run --rm -e wkdir=/data bip   bash /BIP/SCRIPTS/BIP.sh   --platform nanopore --minreads 2 BAI