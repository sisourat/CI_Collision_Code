Steps to run a CI Coll job (no capture)

1) Run Gamess-US to get the Target data: MO, HF energies, integrals for CI
2) Convert Gamess-US output files into Coll input files
   python /home/nicolas/progs/Coll/Coll_with_gus/Python_scripts/from_gus_to_coll_input.py 0

3a) Run CI calculations to get Target states (see MyCIcode)
    As output files one gets f_cistat which containes ncistate CI states

3b) TO BE UPDATED (NOW ONE SHOULD USED GETSTA-INTS-CODE ...)  

    ############Run /home/nicolas/progs/Coll/Coll_with_gus/noprop_sources/Coll to get the V_P matrix elements in the target MO basis set
    ############(same inputs as previous Coll code)
    ############As output files one gets matcoupbI matovlbI where I goes from 1 to nb of impact parameters

4) Run matColl to get the final Coupling matrix
    /home/nicolas/progs/Coll/Coll_with_gus/cicode/matColl < input_dyn
    input_dyn should be **
    matcoupb1 matovlb1
    ncistate ncistate f_cistat f_cistat log
    **
    As output file one gets matCI_b

5) Finally run ../../../prop_sources/Prop matCI_b ista (ista is the index of the initial state, usually 1=target GS)

