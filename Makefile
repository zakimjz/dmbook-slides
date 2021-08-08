all: ychap1.pdf ychap2.pdf ychap3.pdf ychap4.pdf ychap5.pdf ychap6.pdf ychap7.pdf ychap8.pdf ychap9.pdf ychap10.pdf ychap11.pdf ychap12.pdf ychap13.pdf ychap14.pdf ychap15.pdf ychap16.pdf ychap17.pdf ychap18.pdf ychap19.pdf ychap20.pdf ychap21.pdf ychap22.pdf ychap23.pdf ychap24.pdf ychap25.pdf ychap26.pdf ychap27.pdf clean

ychap1.pdf: EDA/data/ydata.tex
	rm -f tmpyslides.*
	echo "\def\doxdata{1} \input{yslides}" > tmpyslides.tex
	latexmk -dvi tmpyslides.tex ;  dvipdf tmpyslides; mv tmpyslides.pdf ychap1.pdf

ychap2.pdf: EDA/numeric/ynumeric.tex
	rm -f tmpyslides.*
	echo "\def\doxnum{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides; mv tmpyslides.pdf ychap2.pdf

ychap3.pdf: EDA/categorical/ycategorical.tex
	rm -f tmpyslides.*
	echo "\def\doxcat{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap3.pdf

ychap4.pdf: EDA/graph/ygraph.tex
	rm -f tmpyslides.*
	echo "\def\doxgraph{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap4.pdf

ychap5.pdf: EDA/kernel/ykernel.tex
	rm -f tmpyslides.*
	echo "\def\doxkernel{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap5.pdf

ychap6.pdf: EDA/highdim/yhighdim.tex
	rm -f tmpyslides.*
	echo "\def\doxhighdim{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap6.pdf

ychap7.pdf: EDA/dimreduction/ydimreduction.tex
	rm -f tmpyslides.*
	echo "\def\doxdimred{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap7.pdf

ychap8.pdf: FPM/itemsets/yitemsets.tex
	rm -f tmpyslides.*
	echo "\def\doxitemsets{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap8.pdf

ychap9.pdf: FPM/sumrep/ysumrep.tex
	rm -f tmpyslides.*
	echo "\def\doxsumrep{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap9.pdf

ychap10.pdf: FPM/sequences/ysequences.tex
	rm -f tmpyslides.*
	echo "\def\doxsequences{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ; dvipdf tmpyslides;   mv tmpyslides.pdf ychap10.pdf

ychap11.pdf: FPM/graphs/ygraphs.tex
	rm -f tmpyslides.*
	echo "\def\doxgraphs{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ; dvipdf tmpyslides;   mv tmpyslides.pdf ychap11.pdf

ychap12.pdf: FPM/fpmeval/yfpmeval.tex
	rm -f tmpyslides.*
	echo "\def\doxfpmeval{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap12.pdf

ychap13.pdf: CLUST/representative/yrepresentative.tex
	rm -f tmpyslides.*
	echo "\def\doxrep{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap13.pdf

ychap14.pdf: CLUST/hierarchical/yhierarchical.tex
	rm -f tmpyslides.*
	echo "\def\doxhier{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap14.pdf

ychap15.pdf: CLUST/density/ydensity.tex
	rm -f tmpyslides.*
	echo "\def\doxdens{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi -latex='latex --shell-escape %O %S' tmpyslides 
	dvips tmpyslides;  ps2pdf -dNOSAFER tmpyslides.ps tmpyslides.pdf
	mv tmpyslides.pdf ychap15.pdf

ychap16.pdf: CLUST/spectral/yspectral.tex
	rm -f tmpyslides.*
	echo "\def\doxspectral{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap16.pdf

ychap17.pdf: CLUST/eval/yeval.tex
	rm -f tmpyslides.*
	echo "\def\doxclusteval{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap17.pdf

ychap18.pdf: CLASS/probabilistic/yprobabilistic.tex
	rm -f tmpyslides.*
	echo "\def\doxprob{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap18.pdf

ychap19.pdf: CLASS/decisiontrees/ydecisiontrees.tex
	rm -f tmpyslides.*
	echo "\def\doxdectrees{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ; dvipdf tmpyslides;   mv tmpyslides.pdf ychap19.pdf

ychap20.pdf: CLASS/lda/ylda.tex
	rm -f tmpyslides.*
	echo "\def\doxlda{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap20.pdf

ychap21.pdf: CLASS/svm/ysvm.tex
	rm -f tmpyslides.*
	echo "\def\doxsvm{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap21.pdf

ychap22.pdf: CLASS/eval/yeval.tex
	rm -f tmpyslides.*
	echo "\def\doxclasseval{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap22.pdf

ychap23.pdf: REG/linear/ylinear.tex
	rm -f tmpyslides.*
	echo "\def\doxlinear{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap23.pdf

ychap24.pdf: REG/logit/ylogit.tex
	rm -f tmpyslides.*
	echo "\def\doxlogit{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap24.pdf

ychap25.pdf: REG/neural/yneural.tex
	rm -f tmpyslides.*
	echo "\def\doxneural{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap25.pdf

ychap26.pdf: REG/deep/ydeep.tex
	rm -f tmpyslides.*
	echo "\def\doxdeep{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap26.pdf

ychap27.pdf: REG/eval/yeval.tex
	rm -f tmpyslides.*
	echo "\def\doxevalreg{1} \input{yslides}">tmpyslides.tex
	latexmk -dvi tmpyslides ;  dvipdf tmpyslides;  mv tmpyslides.pdf ychap27.pdf

clean:
	rm -f tmpyslides.*
	find . -iname "*.aux" -print | xargs rm
