# coding: utf-8
# Input Sequence is the candidate primer_seq and ref_seq, select the probe sequence based on the following principles:
# 1. no G at 5' end of the primer_seq.
# 2. The No. of C is larger than the No. of G in probe.
# 3. The Tm of the probe is 10C higher than the primers. 70~75C for non MGB probe, 60~65C for MGB probe.
# 4. The degenerate number of the probe should be minimal.
# 5. No interactions between the primer and probe, no 3' interaction calculation for probe.
