
REPS= mu1_6 mu2 #mu0 mu0_4 mu0_8 mu1_2

.SILENT: all run tidy clean

all: run tidy

run:
	for r in $(REPS); do \
		echo "$$r"; \
		cd "$$r" && make --silent;\
		cd ..; \
	done

tidy:
	for r in $(REPS); do \
		cd "$$r" && make tidy --silent ;\
		cd ..; \
	done

clean:
	for r in $(REPS); do \
		cd "$$r" && make clean --silent;\
		cd ..; \
	done




.PHONY: run
