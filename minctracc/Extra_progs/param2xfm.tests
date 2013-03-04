#! /bin/sh

check_transform () {
    ./param2xfm -clobber $@ test1.xfm
    cat > test2.xfm
    if ./cmpxfm -linear_tolerance 0.0001 -translation_tolerance 0.0001 \
	test1.xfm test2.xfm; then
	:
    else
	echo >&2 $0 failed: param2xfm $@ produced incorrect results.
	exit 1
    fi
}


check_transform -translation 2 3 4 <<EOF
MNI Transform File

Transform_Type = Linear;
Linear_Transform =
 1 0 0 2
 0 1 0 3
 0 0 1 4;
EOF


check_transform -rotations 0 0 90 <<EOF
MNI Transform File

Transform_Type = Linear;
Linear_Transform =
 0 -1 0 0
 1  0 0 0
 0  0 1 0;
EOF


check_transform -scales 2 3 4 <<EOF
MNI Transform File

Transform_Type = Linear;
Linear_Transform =
 2 0 0 0
 0 3 0 0
 0 0 4 0;
EOF


check_transform -shears 2 3 4 <<EOF
MNI Transform File

Transform_Type = Linear;
Linear_Transform =
 1 0 0 0
 2 1 0 0
 3 4 1 0;
EOF


