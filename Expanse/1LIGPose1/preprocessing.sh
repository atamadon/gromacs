#!/bin/bash

PROTEIN="structures/protein.pdb"
LIGAND_MOL="structures/ligand.mol2"

LIGAND="MMFF/ligand.pdb"
# LIGAND_TOPOLOGY="MMFF/ligand.itp"
LIGAND_TOPOLOGY="ligand.itp"
# LIGAND_PARAMETERS="structures/ligand.prm"
LIGAND_PARAMETERS="ligand.prm"
LIGAND_NAME="LIG"

FORCE_FIELD="charmm36-jul2022.ff"

# CHARMM36 force field
wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz
tar -zxvf *.tgz
rm *.tgz

# ions.mdp
wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp

gmx pdb2gmx -f $PROTEIN -ignh -o conf.pdb << EOF
1
1
EOF

RESPONSE=$(curl -s -F "myMol2=@${LIGAND_MOL}" "swissparam.ch:5678/startparam?approach=mmff-based")
SESSION_NUMBER=$(echo "$RESPONSE" | tr -d '\0' | grep -oP 'sessionNumber=\K[0-9]+')
if [ -z "$SESSION_NUMBER" ]; then
    echo "Failed to retrieve session number. Response: $RESPONSE"
    exit 1
fi
echo "Job submitted. Session number: ${SESSION_NUMBER}"
sleep 10

# Retrieve the results
curl -s "swissparam.ch:5678/retrievesession?sessionNumber=${SESSION_NUMBER}" -o results.tar.gz

if [ $? -eq 0 ]; then
    echo "Results retrieved and saved as results.tar.gz"
else
    echo "Failed to retrieve results."
    exit 1
fi
mkdir MMFF && tar -zxvf results.tar.gz -C MMFF
rm results.tar.gz

awk '/\[ atomtypes \]/{flag=1}/\[ pairtypes \]/{flag=1}/^\[/{if (!/atomtypes|pairtypes/) flag=0}flag' MMFF/ligand.itp > $LIGAND_PARAMETERS

sed -e '/\[ atomtypes \]/,/^\[/{//!d; /\[ atomtypes \]/d}' -e '/\[ pairtypes \]/,/^\[/{//!d; /\[ pairtypes \]/d}' MMFF/ligand.itp > $LIGAND_TOPOLOGY

# Pause the script and wait for user input to perform the following steps

# After generating the topol.top file, insert the following line after the Position restraint file
# #include "PIPC.itp"

# Under the [ molecules ] section, add the following line:
# PIPC 1

# Once the step are complete and the topol.top file is updated, continue with the script
# read -p "Please edit topol.top, press Enter to continue..."

sed -i '/#include ".\/charmm36-jul2022.ff\/forcefield.itp"/a #include "'"$LIGAND_PARAMETERS"'"' topol.top
sed -i '/#include "topol_Protein_chain_B.itp"/a #include "'"$LIGAND_TOPOLOGY"'"' topol.top

echo "${LIGAND_NAME} 1" >> topol.top

# read -p "Please edit topol.top, press Enter to continue..."
cat conf.pdb $LIGAND | grep -v END > system.pdb

gmx editconf -f system.pdb -c -d 1.5 -o boxed.gro 

gmx solvate -cs -cp boxed.gro -p topol.top -o solvated.gro

# read -p "Please edit topol.top, press Enter to continue..."

gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 << EOF
15
EOF

gmx make_ndx -f solv_ions.gro << EOF
1 | 20
q
EOF

echo "System prepared successfully! Proceed with Energy Minimization and Equilibration steps."