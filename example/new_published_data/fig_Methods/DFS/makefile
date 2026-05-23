

MCD_SOLVER      ?= ../../../../build/mcd_solver
QMATRIX_PY      ?= ../../../../build/qMatrixGreenCompute.py
AV_COMBINE_PY   ?= ../../../../build/av_combine.py

SUFF ?= 0.1.txt 0.2.txt #0.3.txt #0.txt 0.05.txt

FILES ?= pFiles/

all:
	mkdir -p logs greens qMatrices; \
	for s in $(SUFF); do \
		echo "Doing $(FILES)pFile_05_$$s"; \
		$(MCD_SOLVER) "$(FILES)pFile_05_$$s" 1 > "logs/log_$$s"; \
		$(QMATRIX_PY) qMatrices.txt -10 10; \
		mv greenFunctionValueQM.txt "green_05_$$s"; \
		$(AV_COMBINE_PY) "green_05_$$s"; \
		mv "green_05_$$s" greens/; \
		mv qMatrices.txt "qMatrices/qM_05_$$s"; \
		mv "av_green_05_$$s" "dos_green$$s"; \
	done; \

tidy:
	rm -fr greens logs qMatrices;

clean:
	rm -fr greens logs qMatrices;
	rm -f dos_green*;
