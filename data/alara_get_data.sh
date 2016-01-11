user_agent='Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.6) Gecko/20070802 SeaMonkey/1.1.4'
rm -f fendlg-2.0_175-gz
wget -c -nc -U "$user_agent" https://www-nds.iaea.org/fendl2/activation/processed/vitj_e/libout/fendlg-2.0_175-gz
rm -f fendld_2.0
touch fendld_2.0
for decay_num in 01 02 03 04 05 06 07 08 09 10; do 
   decay_file=fendld_2.0${decay_num}-gz
   wget -c -nc -U "$user_agent" https://www-nds.iaea.org/fendl2/decay/fendld/${decay_file}
done
