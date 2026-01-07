
REPS= 4s_u8 8s_u8 12s_u8 16s_u8

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
