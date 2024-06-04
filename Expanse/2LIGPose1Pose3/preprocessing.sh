#!/bin/bash

PROTEIN="structures/protein.pdb"

LIGAND_1="../SwissParam/Ligand1/ligand1.pdb"
# LIGAND_2="../SwissParam/Ligand2/ligand2.pdb"
LIGAND_3="../SwissParam/Ligand3/ligand3.pdb"
LIGAND="${LIGAND_1} ${LIGAND_2} ${LIGAND_3}"
NUM_LIGANDS=2

LIGAND_TOPOLOGY_1="../SwissParam/Ligand1/ligand1.itp"
LIGAND_TOPOLOGY_2="../SwissParam/Ligand2/ligand2.itp"
LIGAND_TOPOLOGY_3="../SwissParam/Ligand3/ligand3.itp"

LIGAND_PARAMETERS_1="../SwissParam/Ligand1/ligand1.prm"
LIGAND_PARAMETERS_2="../SwissParam/Ligand2/ligand2.prm"
LIGAND_PARAMETERS_3="../SwissParam/Ligand3/ligand3.prm"

# LIGAND_PARAMETERS="ligand.prm"
LIGAND_NAME="LIG"
FORCE_FIELD="charmm36-jul2022.ff"

# CHARMM36 force field
wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz
tar -zxvf *.tgz
rm *.tgz
# ions.mdp
wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp

# Create topology
gmx pdb2gmx -f $PROTEIN -ignh -o conf.pdb << EOF
1
1
EOF
# Submit ligand to SwissParam
# RESPONSE=$(curl -s -F "myMol2=@${LIGAND_MOL}" "swissparam.ch:5678/startparam?approach=mmff-based")
# SESSION_NUMBER=$(echo "$RESPONSE" | tr -d '\0' | grep -oP 'sessionNumber=\K[0-9]+')
# if [ -z "$SESSION_NUMBER" ]; then
#     echo "Failed to retrieve session number. Response: $RESPONSE"
#     exit 1
# fi
# echo "Job submitted. Session number: ${SESSION_NUMBER}"
# sleep 10
# # Retrieve the results
# curl -s "swissparam.ch:5678/retrievesession?sessionNumber=${SESSION_NUMBER}" -o results.tar.gz
# if [ $? -eq 0 ]; then
#     echo "Results retrieved and saved as results.tar.gz"
# else
#     echo "Failed to retrieve results."
#     exit 1
# fi
# mkdir MMFF && tar -zxvf results.tar.gz -C MMFF
# rm results.tar.gz

# Create ligand parameters and new topology
# awk '/\[ atomtypes \]/{flag=1}/\[ pairtypes \]/{flag=1}/^\[/{if (!/atomtypes|pairtypes/) flag=0}flag' MMFF/ligand.itp > $LIGAND_PARAMETERS
# sed -e '/\[ atomtypes \]/,/^\[/{//!d; /\[ atomtypes \]/d}' -e '/\[ pairtypes \]/,/^\[/{//!d; /\[ pairtypes \]/d}' MMFF/ligand.itp > $LIGAND_TOPOLOGY

# Add ligand topology and parameters to the system topology
# for ((i=1; i<=$NUM_LIGANDS; i++))
# do
#     LIGAND_PARAMETERS="LIGAND_PARAMETERS_$i"
#     LIGAND_TOPOLOGY="LIGAND_TOPOLOGY_$i"
#     sed -i '/#include ".\/charmm36-jul2022.ff\/forcefield.itp"/a #include "'"${!LIGAND_PARAMETERS}"'"' topol.top
#     sed -i '/#include "topol_Protein_chain_B.itp"/a #include "'"${!LIGAND_TOPOLOGY}"'"' topol.top
# done
sed -i '/#include ".\/charmm36-jul2022.ff\/forcefield.itp"/a #include "'"$LIGAND_PARAMETERS_1"'"' topol.top
sed -i '/#include "topol_Protein_chain_B.itp"/a #include "'"$LIGAND_TOPOLOGY_1"'"' topol.top
echo "${LIGAND_NAME} ${NUM_LIGANDS}" >> topol.top

# Create system
cat conf.pdb $LIGAND | grep -v END > system.pdb
# Create box
gmx editconf -f system.pdb -c -d 1.5 -o boxed.gro 
# Solvate
gmx solvate -cs -cp boxed.gro -p topol.top -o solvated.gro
# Add ions
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 << EOF
SOL
EOF
# Prepare system
gmx make_ndx -f solv_ions.gro << EOF
1 | 20
q
EOF

echo "System prepared successfully! Proceed with Energy Minimization and Equilibration steps."