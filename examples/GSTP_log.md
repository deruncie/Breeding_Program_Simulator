Year 1.00
- Created PYT from none.
  - Will be available Year 2.00.
  - Selection=Bootstrap PYT from initial DH seed.
- Genotyped PYT:Year1.00 population using snpChip=1 (n=120).
  - Will be available Year 1.50.

Year 2.00
- Created RC from PYT:Year1.00.
  - Will be available Year 2.00.
  - Selection=Random sample from latest PYT to initialize RC.
- Created DH_PIPE from RC:Year2.00.
  - Will be available Year 4.00.
  - Selection=Random selection from RC; Crossing=makeDH(nDH=1)
    per selected line.
- Started a AYT trial by selecting n=40 from PYT:Year1.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 2.50.

Year 2.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year1.00
  (n_total=120, trait=1, chip=1).
- Created RC from RC:Year2.00.
  - Will be available Year 2.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year2.50 population using snpChip=1 (n=120).
  - Will be available Year 2.75.

Year 2.75
- Created RC from RC:Year2.50.
  - Will be available Year 3.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year2.75 population using snpChip=1 (n=120).
  - Will be available Year 3.00.

Year 3.00
- Created DH_PIPE from RC:Year2.75.
  - Will be available Year 5.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a AYT trial by selecting n=40 from PYT:Year1.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 3.50.
- Started a EYT trial by selecting n=6 from AYT:Year2.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 4.50.
- Created RC from RC:Year2.75.
  - Will be available Year 3.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year3.00 population using snpChip=1 (n=120).
  - Will be available Year 3.25.

Year 3.25
- Created RC from RC:Year3.00.
  - Will be available Year 3.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year3.25 population using snpChip=1 (n=120).
  - Will be available Year 3.50.

Year 3.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year1.00
  (n_total=120, trait=1, chip=1).
- Created RC from RC:Year3.25.
  - Will be available Year 3.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year3.50 population using snpChip=1 (n=120).
  - Will be available Year 3.75.

Year 3.75
- Created RC from RC:Year3.50.
  - Will be available Year 4.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year3.75 population using snpChip=1 (n=120).
  - Will be available Year 4.00.

Year 4.00
- Created DH_PIPE from RC:Year3.75.
  - Will be available Year 6.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year2.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 4.50.
- Genotyped PYT:Year4.00 population using snpChip=1 (n=120).
  - Will be available Year 4.50.
- Started a AYT trial by selecting n=40 from PYT:Year1.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 4.50.
- Started a EYT trial by selecting n=6 from AYT:Year3.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 5.50.
- Created RC from RC:Year3.75.
  - Will be available Year 4.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year4.00 population using snpChip=1 (n=120).
  - Will be available Year 4.25.

Year 4.25
- Created RC from RC:Year4.00.
  - Will be available Year 4.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year4.25 population using snpChip=1 (n=120).
  - Will be available Year 4.50.

Year 4.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year4.00,
  PYT:Year1.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year4.25.
  - Will be available Year 4.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year4.50 population using snpChip=1 (n=120).
  - Will be available Year 4.75.

Year 4.75
- Created RC from RC:Year4.50.
  - Will be available Year 5.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year4.75 population using snpChip=1 (n=120).
  - Will be available Year 5.00.

Year 5.00
- Created DH_PIPE from RC:Year4.75.
  - Will be available Year 7.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year3.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 5.50.
- Genotyped PYT:Year5.00 population using snpChip=1 (n=120).
  - Will be available Year 5.50.
- Started a AYT trial by selecting n=40 from PYT:Year4.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 5.50.
- Started a EYT trial by selecting n=6 from AYT:Year4.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 6.50.
- Created Variety from EYT:Year3.00.
  - Will be available Year 5.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year4.75.
  - Will be available Year 5.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year5.00 population using snpChip=1 (n=120).
  - Will be available Year 5.25.

Year 5.25
- Created RC from RC:Year5.00.
  - Will be available Year 5.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year5.25 population using snpChip=1 (n=120).
  - Will be available Year 5.50.

Year 5.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year5.00,
  PYT:Year4.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year5.25.
  - Will be available Year 5.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year5.50 population using snpChip=1 (n=120).
  - Will be available Year 5.75.

Year 5.75
- Created RC from RC:Year5.50.
  - Will be available Year 6.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year5.75 population using snpChip=1 (n=120).
  - Will be available Year 6.00.

Year 6.00
- Created DH_PIPE from RC:Year5.75.
  - Will be available Year 8.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year4.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 6.50.
- Genotyped PYT:Year6.00 population using snpChip=1 (n=120).
  - Will be available Year 6.50.
- Started a AYT trial by selecting n=40 from PYT:Year5.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 6.50.
- Started a EYT trial by selecting n=6 from AYT:Year5.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 7.50.
- Created Variety from EYT:Year4.00.
  - Will be available Year 6.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year5.75.
  - Will be available Year 6.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year6.00 population using snpChip=1 (n=120).
  - Will be available Year 6.25.

Year 6.25
- Created RC from RC:Year6.00.
  - Will be available Year 6.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year6.25 population using snpChip=1 (n=120).
  - Will be available Year 6.50.

Year 6.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year6.00,
  PYT:Year5.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year6.25.
  - Will be available Year 6.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year6.50 population using snpChip=1 (n=120).
  - Will be available Year 6.75.

Year 6.75
- Created RC from RC:Year6.50.
  - Will be available Year 7.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year6.75 population using snpChip=1 (n=120).
  - Will be available Year 7.00.

Year 7.00
- Created DH_PIPE from RC:Year6.75.
  - Will be available Year 9.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year5.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 7.50.
- Genotyped PYT:Year7.00 population using snpChip=1 (n=120).
  - Will be available Year 7.50.
- Started a AYT trial by selecting n=40 from PYT:Year6.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 7.50.
- Started a EYT trial by selecting n=6 from AYT:Year6.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 8.50.
- Created Variety from EYT:Year5.00.
  - Will be available Year 7.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year6.75.
  - Will be available Year 7.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year7.00 population using snpChip=1 (n=120).
  - Will be available Year 7.25.

Year 7.25
- Created RC from RC:Year7.00.
  - Will be available Year 7.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year7.25 population using snpChip=1 (n=120).
  - Will be available Year 7.50.

Year 7.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year7.00,
  PYT:Year6.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year7.25.
  - Will be available Year 7.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year7.50 population using snpChip=1 (n=120).
  - Will be available Year 7.75.

Year 7.75
- Created RC from RC:Year7.50.
  - Will be available Year 8.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year7.75 population using snpChip=1 (n=120).
  - Will be available Year 8.00.

Year 8.00
- Created DH_PIPE from RC:Year7.75.
  - Will be available Year 10.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year6.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 8.50.
- Genotyped PYT:Year8.00 population using snpChip=1 (n=120).
  - Will be available Year 8.50.
- Started a AYT trial by selecting n=40 from PYT:Year7.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 8.50.
- Started a EYT trial by selecting n=6 from AYT:Year7.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 9.50.
- Created Variety from EYT:Year6.00.
  - Will be available Year 8.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year7.75.
  - Will be available Year 8.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year8.00 population using snpChip=1 (n=120).
  - Will be available Year 8.25.

Year 8.25
- Created RC from RC:Year8.00.
  - Will be available Year 8.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year8.25 population using snpChip=1 (n=120).
  - Will be available Year 8.50.

Year 8.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year8.00,
  PYT:Year7.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year8.25.
  - Will be available Year 8.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year8.50 population using snpChip=1 (n=120).
  - Will be available Year 8.75.

Year 8.75
- Created RC from RC:Year8.50.
  - Will be available Year 9.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year8.75 population using snpChip=1 (n=120).
  - Will be available Year 9.00.

Year 9.00
- Created DH_PIPE from RC:Year8.75.
  - Will be available Year 11.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year7.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 9.50.
- Genotyped PYT:Year9.00 population using snpChip=1 (n=120).
  - Will be available Year 9.50.
- Started a AYT trial by selecting n=40 from PYT:Year8.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 9.50.
- Started a EYT trial by selecting n=6 from AYT:Year8.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 10.50.
- Created Variety from EYT:Year7.00.
  - Will be available Year 9.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year8.75.
  - Will be available Year 9.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year9.00 population using snpChip=1 (n=120).
  - Will be available Year 9.25.

Year 9.25
- Created RC from RC:Year9.00.
  - Will be available Year 9.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year9.25 population using snpChip=1 (n=120).
  - Will be available Year 9.50.

Year 9.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year9.00,
  PYT:Year8.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year9.25.
  - Will be available Year 9.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year9.50 population using snpChip=1 (n=120).
  - Will be available Year 9.75.

Year 9.75
- Created RC from RC:Year9.50.
  - Will be available Year 10.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year9.75 population using snpChip=1 (n=120).
  - Will be available Year 10.00.

Year 10.00
- Created DH_PIPE from RC:Year9.75.
  - Will be available Year 12.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year8.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 10.50.
- Genotyped PYT:Year10.00 population using snpChip=1 (n=120).
  - Will be available Year 10.50.
- Started a AYT trial by selecting n=40 from PYT:Year9.00 by Top
  by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 10.50.
- Started a EYT trial by selecting n=6 from AYT:Year9.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 11.50.
- Created Variety from EYT:Year8.00.
  - Will be available Year 10.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year9.75.
  - Will be available Year 10.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year10.00 population using snpChip=1 (n=120).
  - Will be available Year 10.25.

Year 10.25
- Created RC from RC:Year10.00.
  - Will be available Year 10.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year10.25 population using snpChip=1 (n=120).
  - Will be available Year 10.50.

Year 10.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year10.00,
  PYT:Year9.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year10.25.
  - Will be available Year 10.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year10.50 population using snpChip=1 (n=120).
  - Will be available Year 10.75.

Year 10.75
- Created RC from RC:Year10.50.
  - Will be available Year 11.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year10.75 population using snpChip=1 (n=120).
  - Will be available Year 11.00.

Year 11.00
- Created DH_PIPE from RC:Year10.75.
  - Will be available Year 13.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year9.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 11.50.
- Genotyped PYT:Year11.00 population using snpChip=1 (n=120).
  - Will be available Year 11.50.
- Started a AYT trial by selecting n=40 from PYT:Year10.00 by
  Top by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 11.50.
- Started a EYT trial by selecting n=6 from AYT:Year10.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 12.50.
- Created Variety from EYT:Year9.00.
  - Will be available Year 11.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year10.75.
  - Will be available Year 11.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year11.00 population using snpChip=1 (n=120).
  - Will be available Year 11.25.

Year 11.25
- Created RC from RC:Year11.00.
  - Will be available Year 11.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year11.25 population using snpChip=1 (n=120).
  - Will be available Year 11.50.

Year 11.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year11.00,
  PYT:Year10.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year11.25.
  - Will be available Year 11.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year11.50 population using snpChip=1 (n=120).
  - Will be available Year 11.75.

Year 11.75
- Created RC from RC:Year11.50.
  - Will be available Year 12.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year11.75 population using snpChip=1 (n=120).
  - Will be available Year 12.00.

Year 12.00
- Created DH_PIPE from RC:Year11.75.
  - Will be available Year 14.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year10.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 12.50.
- Genotyped PYT:Year12.00 population using snpChip=1 (n=120).
  - Will be available Year 12.50.
- Started a AYT trial by selecting n=40 from PYT:Year11.00 by
  Top by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 12.50.
- Started a EYT trial by selecting n=6 from AYT:Year11.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 13.50.
- Created Variety from EYT:Year10.00.
  - Will be available Year 12.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year11.75.
  - Will be available Year 12.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year12.00 population using snpChip=1 (n=120).
  - Will be available Year 12.25.

Year 12.25
- Created RC from RC:Year12.00.
  - Will be available Year 12.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year12.25 population using snpChip=1 (n=120).
  - Will be available Year 12.50.

Year 12.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year12.00,
  PYT:Year11.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year12.25.
  - Will be available Year 12.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year12.50 population using snpChip=1 (n=120).
  - Will be available Year 12.75.

Year 12.75
- Created RC from RC:Year12.50.
  - Will be available Year 13.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year12.75 population using snpChip=1 (n=120).
  - Will be available Year 13.00.

Year 13.00
- Created DH_PIPE from RC:Year12.75.
  - Will be available Year 15.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year11.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 13.50.
- Genotyped PYT:Year13.00 population using snpChip=1 (n=120).
  - Will be available Year 13.50.
- Started a AYT trial by selecting n=40 from PYT:Year12.00 by
  Top by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 13.50.
- Started a EYT trial by selecting n=6 from AYT:Year12.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 14.50.
- Created Variety from EYT:Year11.00.
  - Will be available Year 13.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year12.75.
  - Will be available Year 13.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year13.00 population using snpChip=1 (n=120).
  - Will be available Year 13.25.

Year 13.25
- Created RC from RC:Year13.00.
  - Will be available Year 13.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year13.25 population using snpChip=1 (n=120).
  - Will be available Year 13.50.

Year 13.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year13.00,
  PYT:Year12.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year13.25.
  - Will be available Year 13.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year13.50 population using snpChip=1 (n=120).
  - Will be available Year 13.75.

Year 13.75
- Created RC from RC:Year13.50.
  - Will be available Year 14.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year13.75 population using snpChip=1 (n=120).
  - Will be available Year 14.00.

Year 14.00
- Created DH_PIPE from RC:Year13.75.
  - Will be available Year 16.00.
  - Selection=Top by EBV from RC; Crossing=makeDH(nDH=1) per
    selected line.
- Started a PYT trial by selecting n=120 from DH_PIPE:Year12.00
  by Random subset from latest DH_PIPE.
  - The trial has 1 locations with 1 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 14.50.
- Genotyped PYT:Year14.00 population using snpChip=1 (n=120).
  - Will be available Year 14.50.
- Started a AYT trial by selecting n=40 from PYT:Year13.00 by
  Top by EBV from latest PYT (fallback: phenotype).
  - The trial has 4 locations with 2 rep per location and takes
    0.50 years to complete and measures traits 1.
  - Will be available Year 14.50.
- Started a EYT trial by selecting n=6 from AYT:Year13.00 by Top
  by phenotype from latest AYT.
  - The trial has 20 locations with 4 rep per location and takes
    1.50 years to complete and measures traits 1.
  - Will be available Year 15.50.
- Created Variety from EYT:Year12.00.
  - Will be available Year 14.00.
  - Selection=Top by phenotype from latest EYT.
- Created RC from RC:Year13.75.
  - Will be available Year 14.25.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year14.00 population using snpChip=1 (n=120).
  - Will be available Year 14.25.

Year 14.25
- Created RC from RC:Year14.00.
  - Will be available Year 14.50.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year14.25 population using snpChip=1 (n=120).
  - Will be available Year 14.50.

Year 14.50
- Trained GS model AlphaSimR::RRBLUP on PYT:Year14.00,
  PYT:Year13.00 (n_total=240, trait=1, chip=1).
- Created RC from RC:Year14.25.
  - Will be available Year 14.75.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year14.50 population using snpChip=1 (n=120).
  - Will be available Year 14.75.

Year 14.75
- Created RC from RC:Year14.50.
  - Will be available Year 15.00.
  - Selection=Top by EBV from latest RC cohort; Crossing=Random
    mating without selfing among selected RC parents.
- Genotyped RC:Year14.75 population using snpChip=1 (n=120).
  - Will be available Year 15.00.
