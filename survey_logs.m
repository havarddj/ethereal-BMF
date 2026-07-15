
files := ["out_1_0_100.log", "out_11_0_100.log", "out_1_101_200.log", "out_11_101_200.log", "out_11_201_300.log", "out_11_301_400.log", "out_11_401_500.log", "out_11_501_600.log", "out_11_601_650.log", "out_11_651_700.log", "out_11_701_750.log", "out_11_751_800.log", "out_11_801_810.log", "out_11_811_820.log", "out_11_821_830.log", "out_11_831_840.log", "out_11_841_850.log", "out_11_851_860.log", "out_11_861_870.log", "out_11_871_880.log", "out_11_881_890.log", "out_11_891_900.log", "out_11_901_910.log", "out_11_911_920.log", "out_11_921_930.log", "out_11_931_940.log", "out_11_941_950.log", "out_11_951_960.log", "out_11_955_959.log", "out_11_960_960.log", "out_11_961_970.log", "out_11_971_980.log", "out_11_981_990.log", "out_11_991_1000.log", "out_1_201_300.log", "out_1_301_400.log", "out_1_401_500.log", "out_1_501_600.log", "out_1_601_650.log", "out_1_651_700.log", "out_1_701_750.log", "out_1_751_800.log", "out_1_801_820.log", "out_1_821_820.log", "out_1_841_860.log", "out_1_861_880.log", "out_1_881_900.log", "out_1_901_910.log", "out_1_911_920.log", "out_1_921_930.log", "out_1_931_940.log", "out_1_941_950.log", "out_1_951_960.log", "out_1_961_970.log", "out_1_971_980.log", "out_1_981_990.log", "out_1_991_1000.log", "out_2_0_100.log", "out_2_101_200.log", "out_2_201_300.log", "out_2_301_400.log", "out_2_401_500.log", "out_2_501_550.log", "out_2_551_600.log", "out_2_601_650.log", "out_2_651_700.log", "out_2_701_750.log", "out_2_731_740.log", "out_2_741_750.log", "out_2_751_800.log", "out_2_800_800.log", "out_2_801_820.log", "out_2_821_840.log", "out_2_841_860.log", "out_2_861_880.log", "out_2_881_900.log", "out_2_901_910.log", "out_2_911_920.log", "out_2_921_930.log", "out_2_931_940.log", "out_2_941_950.log", "out_2_951_960.log", "out_2_961_970.log", "out_2_971_980.log", "out_2_981_990.log", "out_2_991_1000.log", "out_3_0_500.log", "out_3_501_700.log", "out_3_701_800.log", "out_3_801_900.log", "out_3_901_1000.log", "out_7_0_300.log", "out_7_301_500.log", "out_7_501_600.log", "out_7_601_700.log", "out_7_701_750.log", "out_7_751_800.log", "out_7_801_820.log", "out_7_821_840.log", "out_7_841_860.log", "out_7_851_860.log", "out_7_861_880.log", "out_7_881_900.log", "out_7_901_910.log", "out_7_911_920.log", "out_7_921_930.log", "out_7_931_940.log", "out_7_941_950.log", "out_7_951_960.log", "out_7_961_970.log", "out_7_971_980.log", "out_7_981_990.log", "out_7_991_1000.log"];


forms_written := [];


for d in [1,2,3,7,11] do 
Sprint(d);
forms_d := [];
for f in files do 

i1 := Index(f,"_");
i2 := Index(f[i1+1..#f],"_");


if f[i1+1..i1+i2-1] eq Sprint(d) then 

r := Read("data/logs/" cat f);

lines := Split(r,"\n");

for i in [1..#lines] do  
	// this indicates a class should have been written to file. 
	// we extract the level and characteristic. 
	l := lines[i];
	if l eq "Wrote class to file!" then 

		p1 := lines[i-2];
		p2 := lines[i-1];

		// if it is, this is not the first class written as this level, so we need not recompute
		if p2 ne "Wrote class to file!" then 
			p := p1 cat p2;

			pp := p[75..#p];

			// this occurs because every now and then the characteristic is 
			// too large and the lines aren't split the way the code expects.
			// this is a bit of a hack but it at least runs and prints stuff. 
			try 
				char := StringToInteger(pp[1..Index(pp," ")-1]);
			catch e
				if pp eq "" then
					print "ERROR" cat l;
					Append(~forms_d,"ERROR" cat l);
				else 
					char := StringToInteger(pp[1..Index(pp,"a")-1]);
				end if;
			end try;
			assert IsPrime(char);

			// now we have to track down the level. we go back from the period until we hit a space.
			per_ind := Index(pp,".");
			spc_ind := per_ind;
			if per_ind eq 0 then 
				print "ERROR" cat p;
				Append(~forms_d,"ERROR" cat p);
			else 
				while pp[spc_ind] ne " " and spc_ind ne 0 do 
					spc_ind := spc_ind - 1;
				end while;
				level := pp[spc_ind+1..spc_ind+Index(pp[spc_ind+1..#pp]," ")-1];
				print level cat ";" cat Sprint(char);
				Append(~forms_d,level cat ";" cat Sprint(char));
			end if;
		else
			// in this case, we are noting that another form was written at this level. which needs to be 
			// mentioned, to make sure nothing was left out anywhere. 
			print level cat ";" cat Sprint(char);
			Append(~forms_d,level cat ";" cat Sprint(char));
		end if;
	end if;
end for;
end if;
end for;
print "";
Append(~forms_written,forms_d);
end for;


// once we have all of these, we check against the data files to see 
// if any levels failed to write, got eaten by another process, etc. 

d := 2;
forms_d := forms_written[Index([1,2,3,7,11],d)];

forms_data := Read("data/nonEis_d" cat Sprint(d) cat ".csv");
forms_data := Split(forms_data,"\n");

level_primes := [];
for ff in forms_data[2..#forms_data] do 
	ss := Split(ff,";");
	s := ss[1] cat ";" cat ss[2];
	Append(~level_primes,s);
end for;


for ff in forms_d do
	if not ff in level_primes then 
		print "problem finding " cat ff;
	end if;
end for;





