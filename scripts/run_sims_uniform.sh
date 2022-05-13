#!/bin/bash

OMEGAS=(1 2.077 6)

MH=(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0)

# Declares

FITFILE="HepatitisD.nex.BUSTEDS-MH.fit"
NUMREPLICATES=15

for w in ${OMEGAS[@]}; do
    echo ""
    echo "Simulating under w3: "$w

    for dh in ${MH[@]}; do
        echo "Simulating DH rate: "$dh
	OUTSUFFIX=.o3_"$w"_dh_"$dh"_th_0_th_si_0
	echo "hyphy simulator.bf --fit $FITFILE --output $FITFILE$OUTSUFFIX --busted.test.omega3 $w --busted.test.delta $dh --busted.test.psi 0.0 --busted.text.psi_islands 0.0 --replicates $NUMREPLICATES"
	hyphy simulator.bf --fit $FITFILE --output $FITFILE$OUTSUFFIX --busted.test.omega3 $w --busted.test.delta $dh --busted.test.psi 0.0 --busted.text.psi_islands 0.0 --replicates $NUMREPLICATES
    done

    for th in ${MH[@]}; do
	echo "Simulating TH rate: "$th
	OUTSUFFIX=.o3_"$w"_dh_0_th_"$th"_th_si_0
        echo "hyphy simulator.bf --fit $FITFILE --output $FITFILE$OUTSUFFIX --busted.test.omega3 $w --busted.test.delta 0.0 --busted.test.psi $th --busted.text.psi_islands 0.0 --replicates $NUMREPLICATES"
        hyphy simulator.bf --fit $FITFILE --output $FITFILE$OUTSUFFIX --busted.test.omega3 $w --busted.test.delta 0.0 --busted.test.psi $th --busted.text.psi_islands 0.0 --replicates $NUMREPLICATES
    done

    for thsi in ${MH[@]}; do
	echo "Simulating THSI rate: "$thsi
	OUTSUFFIX=.o3_"$w"_dh_0_th_0_th_si_"$thsi"
	echo "hyphy simulator.bf --fit $FITFILE --output $FITFILE$OUTSUFFIX --busted.test.omega3 $w --busted.test.delta 0.0 --busted.test.psi 0.0 --busted.text.psi_islands $thsi --replicates $NUMREPLICATES"
        hyphy simulator.bf --fit $FITFILE --output $FITFILE$OUTSUFFIX --busted.test.omega3 $w --busted.test.delta 0.0 --busted.test.psi 0.0 --busted.text.psi_islands $thsi --replicates $NUMREPLICATES
    done

done
exit 0
