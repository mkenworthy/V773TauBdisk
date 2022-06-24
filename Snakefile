rule ABCastrom:
    output:
        "src/data/astrom_ABC_chain.h5"
    script:
        "src/scripts/make_astrom_ABC_chain.py"

rule ABorbitbundle:
    output:
        "src/data/v773_tau_AB.hdf5"
    shell:
        "jupyter nbconvert --to notebook --execute src/scripts/make_v773_tau_AB_bundle_notebook.ipynb"
