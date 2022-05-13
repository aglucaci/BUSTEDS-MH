#!/bin/bash

OMEGAS=(1)

MH=(0.05 0.1 0.25 0.5)

# Declares
NUMREPLICATES=100

FITS=(adh.nex.BUSTEDS-MH.fit camelid.nex.BUSTEDS-MH.fit HepatitisD.nex.BUSTEDS-MH.fit)

for FITFILE in ${FITS[@]}; do
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

done
exit 0
