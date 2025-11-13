#!/bin/bash -l 

headerfile=../../header.txt

modify(){
    echo ${1}
    file=$1
    sed -i '/^[[:space:]]*$/d' ${file}
    sed -i '/!>/d' ${file}
    sed -i '/!-/d' ${file}
    sed -i '/!=/d' ${file}

cat > ${file}.tmp << 'endmsg'
!> QMeCha Open release
!>
!> All rights reserved @2025 Matteo Barborini
!> Released under CC-BY-NC-ND license
!>
!> @Author       : Matteo Barborini
!> @Version      : 0.1
!> @Release date : 24.10.2025
!> @Repository   : github.com/QMeCha/QMeCha_open
!>
endmsg
    cat ${file} >> ${file}.tmp 
    mv ${file}.tmp ${file}
}
export -f modify

find . -iname "*.f90" -exec bash -c 'modify "${0}"' {} \;