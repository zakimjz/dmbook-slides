# dmbook-slides
Latex Sources for the slides for Data Mining and Analysis


1. Run "make" to compile all the chapter slides

2. You can also run make ychapN.pdf to compile slides for chap N

3. You will have to increase the tex memory, otherwise the compilation
will fail for chapters 15, 16, and 17.

Recommended parameters:

pool_size=5000000
main_memory=12000000
extra_mem_bot=6000000
font_mem_size=3000000

The specific set of instructions to increase tex memory depend on your tex installation. 

* For MikTex, you can do the following:
 + initexmf --edit-config-file latex 
 + type in the above parameter values
 + initexmf --dump=latex 

* For texlive, you can try the following:
 + locate and edit your texmf.cnf file to include the above parameters
 + run: sudo fmtutil-sys --all


4. Finally, you will have to install the fonts in texmk-local.tar.gz to get the
fonts right. This also depends on your tex installation.

* For MikTex, you can try the following steps:
 + Extract texmf-local.tar.gz and add it to the "Roots"
in "MikTeX Settings (Admin)". After that "Refresh FNDB" and "Update
Formats". 
 + Next, via command line 
   - execute: initexmf --admin --edit-config-file=updmap.cfg
    and then add the following lines
    Map LaTeX_fonts.map
    Map pop.map
    Map Optimademi.map
    Map ptt.map
   - execute initexmf -u, followed by initexmf --mkmaps 

* For texlive, expand/extract the texmf-local.tar.gz file into your
texmf-local directory. After that you should run the following
commands:
 + sudo texhash
 + sudo updmap-sys --enable Map LaTeX_fonts.map
 + sudo updmap-sys --enable Map pop.map
 + sudo updmap-sys --enable Map Optimademi.map
 + sudo updmap-sys --enable Map ptt.map
 + sudo texhash

