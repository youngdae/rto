#!/bin/bash
#                  ramp scale
# case9              0.0004
# case30             0.0005
# case118            0.002
# case300            0.002
# case1354pegase     0.0025
# case2383wp         0.005
# case9241pegase     0.01

function usage
{
    echo "Usage: gensol.sh [-w str] [-p int] [-m int] [-o int] [-q int] [-l str] [-s int]"
    echo "       -w str: cold, shift_copy, or shift_phase1"
    echo "       -p int: 0 for no phase 1, and 1 for phase 1"
    echo "       -m int: number of maximum iterations"
    echo "       -o int: option file number to use for approximation"
    echo "       -q int: 0 for no QP, 1 for QP"
    echo "       -l str: mumps, ma27, or ma57"
    echo "       -s int: 0 for not loading a solution, and 1 for loading a solution"
}

LO=2
T=(10)
H=(1)
SC=("halfhour_30")
nlp="ipopt"
cutline="0"
cutgen="1"
perturb="0"
pfsolve="0"
C=("9")
RS=("0.01")
LS=("1.0")
QP=0
W="cold"
PHASE1=0
APPROX="no"
ML=1
OPT=2       # option file number for QP solve
LA="mumps"  # linear algebra engine
SOL=0       # whether to load a solution

# Read the options.
options=`getopt w:p:m:o:q:l:s:b:f: $*`
if [ $? -ne 0 ]; then
    usage
    exit 1
fi

set -- $options
for i
do
    case "$i" in
	-w)
            W=$2
	    shift
            shift
            ;;
	-p)
            PHASE1=$2
	    shift
            shift
	    [[ ! $PHASE1 =~ 0|1 ]] && {
	        usage
	        exit 1
	    }
	    ;;
	-m)
	    APPROX="yes"
            ML=$2
	    shift
            shift
	    ;;
	-o)
	    OPT=$2
	    shift
            shift
	    [[ $OPT -le 2 ]] && {
	        usage
		exit 1
	    }
	    ;;
	-q)
	    QP=$2
	    shift
            shift
	    [[ ! $QP =~ 0|1 ]] && {
		usage
		exit 1
	    }
	    ;;
	-l)
	    LA=$2
	    shift
	    shift
	    ;;
    -s)
        SOL=$2
        shift
        shift
        ;;
	-b)
	    perturb=$2
	    shift
	    shift
	    ;;
	-f)
	    pfsolve=$2
	    shift
	    shift
	    ;;
	--)
	    shift
	    break;;
    esac
done

WSTR="cold"
PSTR="phase0"
QSTR="no_qp"

if [ ${W} != "cold" ]; then
    WSTR="warm"
fi

if [ ${PHASE1} -eq 1 ]; then
    PSTR="phase1"
fi

if [ ${QP} -eq 1 ]; then
    QSTR="use_qp"
fi

#if [ "${APPROX}" == "yes" ]; then
#    cp -f ipopt.op2 ipopt.op${OPT}
#    echo "max_iter ${ML}" >> ipopt.op${OPT}
#fi

for c in ${!C[@]}
do
    for t in ${T[@]}
    do
	    for h in ${H[@]}
    	do
	        for sc in ${SC[@]}
    	    do
	    	    for ls in ${LS[@]}
    		    do
                    if [ ${APPROX} != "yes" ]; then
                        if [ ${OPT} == 7 ]; then
                            result="${nlp}_case${C[$c]}_t${t}_h${h}_${sc}_ls${ls}_rs${RS[$c]}_cl${cutline}_cg${cutgen}_pt${perturb}_pf${pfsolve}_${WSTR}_${PSTR}_${QSTR}_${LA}_lsdual"
#                            result="${nlp}_case${C[$c]}_t${t}_h${h}_${sc}_ls${ls}_rs${RS[$c]}_cl${cutline}_cg${cutgen}_pt${perturb}_${WSTR}_${PSTR}_${QSTR}_${LA}_lsdual"
                        else
                            result="${nlp}_case${C[$c]}_t${t}_h${h}_${sc}_ls${ls}_rs${RS[$c]}_cl${cutline}_cg${cutgen}_pt${perturb}_pf${pfsolve}_${WSTR}_${PSTR}_${QSTR}_${LA}"
#                            result="${nlp}_case${C[$c]}_t${t}_h${h}_${sc}_ls${ls}_rs${RS[$c]}_cl${cutline}_cg${cutgen}_pt${perturb}_${WSTR}_${PSTR}_${QSTR}_${LA}"
                        fi
                    else
                        result="${nlp}_case${C[$c]}_t${t}_h${h}_${sc}_ls${ls}_rs${RS[$c]}_cl${cutline}_cg${cutgen}_pt${perturb}_pf${pfsolve}_${WSTR}_${PSTR}_${QSTR}_${LA}_mlim_${ML}"
#                        result="${nlp}_case${C[$c]}_t${t}_h${h}_${sc}_ls${ls}_rs${RS[$c]}_cl${cutline}_cg${cutgen}_pt${perturb}_${WSTR}_${PSTR}_${QSTR}_${LA}_mlim_${ML}"
                    fi

                    rm -f ${result}.txt

                    julia --project=. src/mpc.jl "data/case${C[$c]}" "data/case${C[$c]}/$sc" ${t} ${h} "${ls}" "${RS[$c]}" "${W}" ${OPT} "${cutline}" "${cutgen}" "${perturb}" "${QP}" "${SOL}" "${pfsolve}" "${result}" > ${result}.txt 2>&1 &
		done
	    done
	done
    done
done


