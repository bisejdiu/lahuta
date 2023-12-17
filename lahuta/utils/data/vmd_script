proc draw_interaction {atom1_index atom2_index sphere_resolution} {

    set atom1 [atomselect top "index $atom1_index"]
    set atom2 [atomselect top "index $atom2_index"]

    set coord [$atom1 get {x y z}]
    set coord [string trim $coord "{}"]
    lassign [split $coord] x1 y1 z1

    set coord [$atom2 get {x y z}]
    set coord [string trim $coord "{}"]
    lassign [split $coord] x2 y2 z2

    graphics top cylinder [list $x1 $y1 $z1] [list $x2 $y2 $z2] radius 0.1
}

proc draw_interactions {index_pairs sphere_resolution} {
    foreach pair $index_pairs {
        draw_interaction top [lindex $pair 0] [lindex $pair 1] $sphere_resolution
    }
}
