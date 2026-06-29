SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# aKG
for round in 18 19 21
do
	bash ${SCRIPT_DIR}/bash_per.sh aKG ${round}
done


# nKG
for round in 19 20 22
do
	bash ${SCRIPT_DIR}/bash_per.sh nKG ${round}
done
