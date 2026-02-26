Year 1.00
- Created CROSS_BLOCK from none.
  - Will be available Year 1.00.
  - Selection=Initialize crossing block from founder GV.
- Created CROSS_BLOCK from CROSS_BLOCK:Year1.00.
  - Will be available Year 1.25.
  - Selection=Fallback to incumbent crossing block only.
- Created CROSS_SEED from CROSS_BLOCK:Year1.00.
  - Will be available Year 1.25.
  - Crossing=Random pair sampling without replacement from current CROSS_BLOCK.

Year 1.50
- Created DH_BULK from CROSS_SEED:Year1.00.
  - Will be available Year 2.75.
  - Crossing=makeDH(nDH=100) per F1 family.

Year 2.00
- Created CROSS_BLOCK from CROSS_BLOCK:Year1.00.
  - Will be available Year 2.25.
  - Selection=Fallback to incumbent crossing block only.
- Created CROSS_SEED from CROSS_BLOCK:Year2.00.
  - Will be available Year 2.25.
  - Crossing=Random pair sampling without replacement from current CROSS_BLOCK.

Year 2.50
- Created DH_BULK from CROSS_SEED:Year2.00.
  - Will be available Year 3.75.
  - Crossing=makeDH(nDH=100) per F1 family.

Year 3.00
- Created CROSS_BLOCK from CROSS_BLOCK:Year2.00.
  - Will be available Year 3.25.
  - Selection=Fallback to incumbent crossing block only.
- Created CROSS_SEED from CROSS_BLOCK:Year3.00.
  - Will be available Year 3.25.
  - Crossing=Random pair sampling without replacement from current CROSS_BLOCK.
- Started a HEADROW_SEL trial by selecting n=500 from DH_BULK:Year1.50 by
  Visual/headrow selection on phenotype.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 4.00.

Year 3.50
- Created DH_BULK from CROSS_SEED:Year3.00.
  - Will be available Year 4.75.
  - Crossing=makeDH(nDH=100) per F1 family.

Year 4.00
- Created CROSS_BLOCK from CROSS_BLOCK:Year3.00.
  - Will be available Year 4.25.
  - Selection=Fallback to incumbent crossing block only.
- Created CROSS_SEED from CROSS_BLOCK:Year4.00.
  - Will be available Year 4.25.
  - Crossing=Random pair sampling without replacement from current CROSS_BLOCK.
- Started a HEADROW_SEL trial by selecting n=500 from DH_BULK:Year2.50 by
  Visual/headrow selection on phenotype.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 5.00.
- Started a PYT trial by selecting n=500 from HEADROW_SEL:Year3.00 by Random
  subset from HEADROW_SEL.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 5.00.

Year 4.50
- Created DH_BULK from CROSS_SEED:Year4.00.
  - Will be available Year 5.75.
  - Crossing=makeDH(nDH=100) per F1 family.

Year 5.00
- Created CROSS_BLOCK from CROSS_BLOCK:Year4.00, PYT:Year4.00.
  - Will be available Year 5.25.
  - Selection=Fallback to incumbent crossing block only.
- Created CROSS_SEED from CROSS_BLOCK:Year5.00.
  - Will be available Year 5.25.
  - Crossing=Random pair sampling without replacement from current CROSS_BLOCK.
- Started a HEADROW_SEL trial by selecting n=500 from DH_BULK:Year3.50 by
  Visual/headrow selection on phenotype.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 6.00.
- Started a PYT trial by selecting n=500 from HEADROW_SEL:Year4.00 by Random
  subset from HEADROW_SEL.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 6.00.
- Started a AYT trial by selecting n=50 from PYT:Year4.00 by Top by phenotype
  from latest PYT.
  - The trial has 4 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 6.00.

Year 5.50
- Created DH_BULK from CROSS_SEED:Year5.00.
  - Will be available Year 6.75.
  - Crossing=makeDH(nDH=100) per F1 family.

Year 6.00
- Created CROSS_BLOCK from CROSS_BLOCK:Year5.00, PYT:Year5.00, AYT:Year5.00.
  - Will be available Year 6.25.
  - Selection=20 PYT + 10 AYT + 20 incumbent non-PYT (EYT-informed when
    available).
- Created CROSS_SEED from CROSS_BLOCK:Year6.00.
  - Will be available Year 6.25.
  - Crossing=Random pair sampling without replacement from current CROSS_BLOCK.
- Started a HEADROW_SEL trial by selecting n=500 from DH_BULK:Year4.50 by
  Visual/headrow selection on phenotype.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 7.00.
- Started a PYT trial by selecting n=500 from HEADROW_SEL:Year5.00 by Random
  subset from HEADROW_SEL.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 7.00.
- Started a AYT trial by selecting n=50 from PYT:Year5.00 by Top by phenotype
  from latest PYT.
  - The trial has 4 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 7.00.
- Started a EYT1 trial by selecting n=10 from AYT:Year5.00 by Top by phenotype
  from latest AYT.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 7.00.

Year 6.50
- Created DH_BULK from CROSS_SEED:Year6.00.
  - Will be available Year 7.75.
  - Crossing=makeDH(nDH=100) per F1 family.

Year 7.00
- Created CROSS_BLOCK from CROSS_BLOCK:Year6.00, PYT:Year6.00, AYT:Year6.00.
  - Will be available Year 7.25.
  - Selection=20 PYT + 10 AYT + 20 incumbent non-PYT (EYT-informed when
    available).
- Created CROSS_SEED from CROSS_BLOCK:Year7.00.
  - Will be available Year 7.25.
  - Crossing=Random pair sampling without replacement from current CROSS_BLOCK.
- Started a HEADROW_SEL trial by selecting n=500 from DH_BULK:Year5.50 by
  Visual/headrow selection on phenotype.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 8.00.
- Started a PYT trial by selecting n=500 from HEADROW_SEL:Year6.00 by Random
  subset from HEADROW_SEL.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 8.00.
- Started a AYT trial by selecting n=50 from PYT:Year6.00 by Top by phenotype
  from latest PYT.
  - The trial has 4 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 8.00.
- Started a EYT1 trial by selecting n=10 from AYT:Year6.00 by Top by phenotype
  from latest AYT.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 8.00.
- Started a EYT2 trial by selecting n=10 from EYT1:Year6.00 by Re-evaluate same
  lines from EYT1.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 8.00.

Year 7.50
- Created DH_BULK from CROSS_SEED:Year7.00.
  - Will be available Year 8.75.
  - Crossing=makeDH(nDH=100) per F1 family.

Year 8.00
- Created CROSS_BLOCK from CROSS_BLOCK:Year7.00, PYT:Year7.00, AYT:Year7.00.
  - Will be available Year 8.25.
  - Selection=20 PYT + 10 AYT + 20 incumbent non-PYT (EYT-informed when
    available).
- Created CROSS_SEED from CROSS_BLOCK:Year8.00.
  - Will be available Year 8.25.
  - Crossing=Random pair sampling without replacement from current CROSS_BLOCK.
- Started a HEADROW_SEL trial by selecting n=500 from DH_BULK:Year6.50 by
  Visual/headrow selection on phenotype.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 9.00.
- Started a PYT trial by selecting n=500 from HEADROW_SEL:Year7.00 by Random
  subset from HEADROW_SEL.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 9.00.
- Started a AYT trial by selecting n=50 from PYT:Year7.00 by Top by phenotype
  from latest PYT.
  - The trial has 4 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 9.00.
- Started a EYT1 trial by selecting n=10 from AYT:Year7.00 by Top by phenotype
  from latest AYT.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 9.00.
- Started a EYT2 trial by selecting n=10 from EYT1:Year7.00 by Re-evaluate same
  lines from EYT1.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 9.00.
- Created Variety from EYT2:Year7.00.
  - Will be available Year 8.00.
  - Selection=Best 2-year mean phenotype from EYT1+EYT2.

Year 8.50
- Created DH_BULK from CROSS_SEED:Year8.00.
  - Will be available Year 9.75.
  - Crossing=makeDH(nDH=100) per F1 family.
Year 8-20 (same as Year 8)

Year 21.00
- Created PI_CAND from PYT:Year20.00, CROSS_BLOCK:Year20.00.
  - Will be available Year 21.00.
  - Selection=Warm-start random sample n=550 from PYT + CROSS_BLOCK.
- Genotyped PYT:Year4.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year5.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year6.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year7.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year8.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year9.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year10.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year11.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year12.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year13.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year14.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year15.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year16.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year17.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year18.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year19.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Genotyped PYT:Year20.00 population using snpChip=1 (n=500).
  - Will be available Year 21.00.
- Trained GS model AlphaSimR::RRBLUP on PYT:Year17.00, PYT:Year18.00,
  PYT:Year19.00, PYT:Year20.00 (n_total=2000, trait=1, chip=1).
- Genotyped PI_CAND:Year21.00 population using snpChip=1 (n=550).
  - Will be available Year 21.00.
- Created PI_CAND from PI_CAND:Year21.00.
  - Will be available Year 21.50.
  - Selection=year_1_cycle_A: split male/female then GS top-100 per sex;
    Crossing=Open-pollination approximation: one random selected male per
    selected female.
- Created PD_DH_INPUT from PI_CAND:Year21.00.
  - Will be available Year 21.50.
  - Selection=year_1_cycle_A: reserve 10 seeds/family for DH.
- Started a HEADROW_SEL trial by selecting n=500 from DH_BULK:Year19.50 by
  Advance fixed n from DH_BULK to headrows.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 22.00.
- Started a PYT trial by selecting n=500 from HEADROW_SEL:Year20.00 by Random
  subset from HEADROW_SEL.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 22.00.
- Started a AYT trial by selecting n=50 from PYT:Year20.00 by Top phenotype
  from latest PYT.
  - The trial has 4 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 22.00.
- Started a EYT1 trial by selecting n=10 from AYT:Year20.00 by Top phenotype
  from latest AYT.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 22.00.
- Started a EYT2 trial by selecting n=10 from EYT1:Year20.00 by Re-evaluate
  EYT1 lines.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 22.00.
- Created Variety from EYT2:Year20.00.
  - Will be available Year 21.00.
  - Selection=Best 2-year EYT mean phenotype.
- Genotyped HEADROW_SEL:Year21.00 population using snpChip=1 (n=500).
  - Will be available Year 21.50.
- Genotyped PYT:Year21.00 population using snpChip=1 (n=500).
  - Will be available Year 21.50.
- Trained GS model AlphaSimR::RRBLUP on PYT:Year17.00, PYT:Year18.00,
  PYT:Year19.00, PYT:Year20.00 (n_total=2000, trait=1, chip=1).

Year 21.50
- Genotyped PI_CAND:Year21.00 population using snpChip=1 (n=3000).
  - Will be available Year 21.50.
- Created PI_CAND from PI_CAND:Year21.00.
  - Will be available Year 22.00.
  - Selection=year_1_cycle_B: split male/female then GS top-100 per sex;
    Crossing=Open-pollination approximation: one random selected male per
    selected female.
- Created PD_DH_INPUT from PI_CAND:Year21.00.
  - Will be available Year 22.00.
  - Selection=year_1_cycle_B: reserve 10 seeds/family for DH.

Year 22.00
- Created DH_BULK from PD_DH_INPUT:Year21.00, PD_DH_INPUT:Year21.50.
  - Will be available Year 23.00.
  - Crossing=makeDH(nDH=31) from one representative seed/family.
- Genotyped PI_CAND:Year21.50 population using snpChip=1 (n=3000).
  - Will be available Year 22.00.
- Created PI_CAND from PI_CAND:Year21.50.
  - Will be available Year 22.50.
  - Selection=year_2_cycle_A: split male/female then GS top-100 per sex;
    Crossing=Open-pollination approximation: one random selected male per
    selected female.
- Created PD_DH_INPUT from PI_CAND:Year21.50.
  - Will be available Year 22.50.
  - Selection=year_2_cycle_A: reserve 10 seeds/family for DH.
- Started a HEADROW_SEL trial by selecting n=500 from DH_BULK:Year20.50 by
  Advance fixed n from DH_BULK to headrows.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 23.00.
- Started a PYT trial by selecting n=500 from HEADROW_SEL:Year21.00 by Random
  subset from HEADROW_SEL.
  - The trial has 1 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 23.00.
- Started a AYT trial by selecting n=50 from PYT:Year21.00 by Top phenotype
  from latest PYT.
  - The trial has 4 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 23.00.
- Started a EYT1 trial by selecting n=10 from AYT:Year21.00 by Top phenotype
  from latest AYT.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 23.00.
- Started a EYT2 trial by selecting n=10 from EYT1:Year21.00 by Re-evaluate
  EYT1 lines.
  - The trial has 8 locations with 1 rep per location and takes 1.00 years to
    complete and measures traits 1.
  - Will be available Year 23.00.
- Created Variety from EYT2:Year21.00.
  - Will be available Year 22.00.
  - Selection=Best 2-year EYT mean phenotype.
- Genotyped HEADROW_SEL:Year22.00 population using snpChip=1 (n=500).
  - Will be available Year 22.50.
- Trained GS model AlphaSimR::RRBLUP on PYT:Year18.00, PYT:Year19.00,
  PYT:Year20.00, PYT:Year21.00 (n_total=2000, trait=1, chip=1).

Year 22.50
- Genotyped PI_CAND:Year22.00 population using snpChip=1 (n=3000).
  - Will be available Year 22.50.
- Created PI_CAND from PI_CAND:Year22.00.
  - Will be available Year 23.00.
  - Selection=year_2_cycle_B: split male/female then GS top-100 per sex;
    Crossing=Open-pollination approximation: one random selected male per
    selected female.
- Created PD_DH_INPUT from PI_CAND:Year22.00.
  - Will be available Year 23.00.
  - Selection=year_2_cycle_B: reserve 10 seeds/family for DH.
Year 22-32 (same as Year 22)
