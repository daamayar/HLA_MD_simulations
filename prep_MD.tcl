## USAGE: vmd -dispdev text -e prep_MD.tcl -args <file_name> <namd_dir>

set pdb_file [lindex $argv 0]
set namd_dir [lindex $argv 1]

package require psfgen
resetpsf

topology ${namd_dir}/top_all36_prot.rtf

set mol [ molecule new ${pdb_file}.pdb ]
set sel [atomselect $mol protein]
$sel moveby [vecinvert [measure center $sel weight mass]]
set chains [lsort -unique [$sel get chain]]

foreach chain $chains {
    puts "Adding protein chain $chain to psfgen"
    set seg ${chain}PRO
    set sel [atomselect $mol "protein and chain $chain"]
    $sel set segid $seg
    $sel writepdb tmp.pdb
    segment $seg { pdb tmp.pdb }
    coordpdb tmp.pdb
}

guesscoord

writepsf ${pdb_file}_autopsf.psf
writepdb ${pdb_file}_autopsf.pdb

package require solvate
solvate ${pdb_file}_autopsf.psf ${pdb_file}_autopsf.pdb -t 20 -o ${pdb_file}_solvate

package require autoionize
autoionize -psf ${pdb_file}_solvate.psf -pdb ${pdb_file}_solvate.pdb -neutralize -o ${pdb_file}_ionized

exit
