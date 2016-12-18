all: ychap1.pdf ychap2.pdf ychap3.pdf ychap4.pdf ychap5.pdf ychap6.pdf ychap7.pdf ychap8.pdf ychap9.pdf ychap10.pdf ychap11.pdf ychap12.pdf ychap13.pdf ychap14.pdf ychap15.pdf ychap16.pdf ychap17.pdf ychap18.pdf ychap19.pdf ychap20.pdf ychap21.pdf ychap22.pdf 

ychap1.pdf: EDA/data/ydata.tex
	latex  "\def\doxdata{1} \input{yslides.tex}" ; latex  "\def\doxdata{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap1.pdf
ychap2.pdf: EDA/numeric/ynumeric.tex
	latex  "\def\doxnum{1} \input{yslides.tex}" ; latex  "\def\doxnum{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap2.pdf
ychap3.pdf: EDA/categorical/ycategorical.tex
	latex  "\def\doxcat{1} \input{yslides.tex}" ; latex  "\def\doxcat{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap3.pdf
ychap4.pdf: EDA/graph/ygraph.tex
	latex  "\def\doxgraph{1} \input{yslides.tex}" ; latex  "\def\doxgraph{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap4.pdf
ychap5.pdf: EDA/kernel/ykernel.tex
	latex  "\def\doxkernel{1} \input{yslides.tex}" ; latex  "\def\doxkernel{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap5.pdf
ychap6.pdf: EDA/highdim/yhighdim.tex
	latex  "\def\doxhighdim{1} \input{yslides.tex}" ; latex  "\def\doxhighdim{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap6.pdf
ychap7.pdf: EDA/dimreduction/ydimreduction.tex
	latex  "\def\doxdimred{1} \input{yslides.tex}" ; latex  "\def\doxdimred{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap7.pdf

ychap8.pdf: FPM/itemsets/yitemsets.tex
	latex  "\def\doxitemsets{1} \input{yslides.tex}" ; latex  "\def\doxitemsets{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap8.pdf
ychap9.pdf: FPM/sumrep/ysumrep.tex
	latex  "\def\doxsumrep{1} \input{yslides.tex}" ; latex  "\def\doxsumrep{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap9.pdf
ychap10.pdf: FPM/sequences/ysequences.tex
	latex  "\def\doxsequences{1} \input{yslides.tex}" ; latex  "\def\doxsequences{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap10.pdf
ychap11.pdf: FPM/graphs/ygraphs.tex
	latex  "\def\doxgraphs{1} \input{yslides.tex}" ; latex  "\def\doxgraphs{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap11.pdf
ychap12.pdf: FPM/fpmeval/yfpmeval.tex
	latex  "\def\doxfpmeval{1} \input{yslides.tex}" ; latex  "\def\doxfpmeval{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap12.pdf

ychap13.pdf: CLUST/representative/yrepresentative.tex
	latex  "\def\doxrep{1} \input{yslides.tex}" ; latex  "\def\doxrep{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap13.pdf
ychap14.pdf: CLUST/hierarchical/yhierarchical.tex
	latex  "\def\doxhier{1} \input{yslides.tex}" ; latex  "\def\doxhier{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap14.pdf
ychap15.pdf: CLUST/density/ydensity.tex
	latex  "\def\doxdens{1} \input{yslides.tex}" ; latex  "\def\doxdens{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap15.pdf
ychap16.pdf: CLUST/spectral/yspectral.tex
	latex  "\def\doxspectral{1} \input{yslides.tex}" ; latex  "\def\doxspectral{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap16.pdf
ychap17.pdf: CLUST/eval/yeval.tex
	latex  "\def\doxclusteval{1} \input{yslides.tex}" ; latex  "\def\doxclusteval{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap17.pdf

ychap18.pdf: CLASS/probabilistic/yprobabilistic.tex
	latex  "\def\doxprob{1} \input{yslides.tex}" ; latex  "\def\doxprob{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap18.pdf
ychap19.pdf: CLASS/decisiontrees/ydecisiontrees.tex
	latex  "\def\doxdectrees{1} \input{yslides.tex}" ; latex  "\def\doxdectrees{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap19.pdf
ychap20.pdf: CLASS/lda/ylda.tex
	latex  "\def\doxlda{1} \input{yslides.tex}" ; latex  "\def\doxlda{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap20.pdf
ychap21.pdf: CLASS/svm/ysvm.tex
	latex  "\def\doxsvm{1} \input{yslides.tex}" ; latex  "\def\doxsvm{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap21.pdf
ychap22.pdf: CLASS/eval/yeval.tex
	latex  "\def\doxclasseval{1} \input{yslides.tex}" ; latex  "\def\doxclasseval{1} \input{yslides.tex}" ;  dvips yslides ; ps2pdf yslides.ps ;  mv yslides.pdf ychap22.pdf
