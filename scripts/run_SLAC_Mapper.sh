SLACMAPPER=/mnt/e/BUSTEDS-MH/scripts/slac-mapper.bf

DATADIR=/mnt/e/BUSTEDS-MH/analysis/13-datasets/SLAC

for file in "$DATADIR"/*.json; do
    output="$file".subs
    echo hyphy $SLACMAPPER $file $output
    hyphy $SLACMAPPER $file $output
done

exit 0
