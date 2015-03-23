#!/bin/sh
cd ../data
grep -i 'DNA binding'  original.entry.list > allTFs.list
grep -i 'DNA-binding'  original.entry.list >> allTFs.list
grep -i 'Transcription regulator'  original.entry.list >> allTFs.list
grep -i 'Transcription factor'  original.entry.list >> allTFs.list
sort allTFs.list | uniq -c | sort -nr  > allTFs.uc.list
sort allTFs.list | uniq > allTFs.u.list

#combined with Trire to get a manually annotated list
sort original.trire_regs | uniq | tail --lines=+3  > trireTFs.list
#cut -f1 -d ' ' allTFs.u.list > allTFs.u2.list
grep -f trireTFs.list original.entry.list > trireTFsF.list
cat trireTFsF.list allTFs.u.list | sort | uniq > AllTrireTFs.list
cat trireTFsF.list allTFs.u.list | sort | uniq -c | sort -nr > AllTrireTFsWC.list
cat trireTFsF.list allTFs.u.list | sort | uniq -d | sort | sort | uniq > AllTrireTFsCommon.list
#and some final outputs to see differences
grep -v -f trireTFsF.list allTFs.u.list > NotInTRIRE.list
grep -v -f allTFs.u.list trireTFsF.list > OnlyInTRIRE.list

#combine with Arvas07 all IPRs to get possible fungal domains
#fromdos original.AF04_acc2ipr.txt  needed to be done before
cut -f2 original.AF04_acc2ipr.txt | sort | uniq > AF04u.list
grep -f AF04u.list allTFs.u.list > AllFungiTFs.list

#AllFungiTFs.list should include everything that is in AllTrireTFsCommon.list
#and as it does it is the final list.

ln -s AllFungiTFs.list interpro.input


cd ../bin