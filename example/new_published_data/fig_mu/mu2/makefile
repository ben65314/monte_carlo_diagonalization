

MCD_SOLVER      ?= ../../../../build/mcd_solver
QMATRIX_PY      ?= ../../../../build/qMatrixGreenCompute.py
AV_COMBINE_PY   ?= ../../../../build/av_combine.py

ifeq ($(wildcard dos_green100.txt),)
	#SUFF ?= Ncst.txt 100.txt
	SUFF ?= 100.txt
else
	#SUFF ?= Ncst.txt 
	SUFF ?= 100.txt 
endif

FILES ?= pFiles/

all:
	mkdir -p logs greens qMatrices; \
	for s in $(SUFF); do \
		echo "Doing $(FILES)pFile_$$s"; \
		$(MCD_SOLVER) "$(FILES)pFile_$$s" 1 > "logs/log_$$s"; \
		$(QMATRIX_PY) qMatrices.txt -10 10; \
		mv greenFunctionValueQM.txt "green_$$s"; \
		$(AV_COMBINE_PY) "green_$$s"; \
		mv "green_$$s" greens/; \
		mv qMatrices.txt "qMatrices/qM_$$s"; \
		mv "av_green_$$s" "dos_green$$s"; \
	done; \

tidy:
	rm -fr greens logs qMatrices;

clean:
	rm -fr greens logs qMatrices;
	rm -f dos_green*;
