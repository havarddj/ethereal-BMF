
AttachSpec("../spec");


function VerifiedData(d)
	K := QuadFld(d);
	ZK := MaximalOrder(K);

	input := "data/verified_d" cat Sprint(d) cat ".csv";
	lines := Split(Read(input), "\n");

	p_le20 := [];
	for l in lines do 
		ll := Split(l,";");
		p := StringToInteger(ll[2]);
		if p le 20 then 
			Append(~p_le20,l);
		end if;
	end for;
	return #lines, #p_le20;
end function;


function AllLifts(d)
	K := QuadFld(d);
	ZK := MaximalOrder(K);

	lmfdb_lifts := Split(Read("sage/d" cat Sprint(d) cat "_liftable.csv"));
	irrat_data := Split(Read("data/irrat_lifts_d" cat Sprint(d) cat "_v2.csv"));

	lmfdb_lifts_ple20 := [];
	irrat_lifts := [];
	irrat_lifts_ple20 := [];

	// gather the 2 < p < 20 lmfdb lifts 
	for l in lmfdb_lifts do 
		ll := Split(l,";");
		if ll[1] ne "level" then 
			p := StringToInteger(ll[2]);
			if p le 20 then 
				Append(~lmfdb_lifts_ple20,l);
			end if;
		end if;
	end for;

	// and now the irrational lifts 
	for l in irrat_data do 
		ll := Split(l,";");
		if ll[1] ne "level" then 
			if ll[3] ne "None found" then 
				p := StringToInteger(ll[2]);
				if p le 20 then 
					Append(~irrat_lifts_ple20,l);
				end if;
				Append(~irrat_lifts,l);
			end if;
		end if;
	end for;

	return #lmfdb_lifts + #irrat_lifts, #lmfdb_lifts_ple20 + #irrat_lifts_ple20;
end function;


function EtherealData(d)
	K := QuadFld(d);
	ZK := MaximalOrder(K);

	eth := Split(Read("data/nonEis_d" cat Sprint(d) cat "_v2.csv"));

	ethereal := [];
	ethereal_ple20 := [];

	for l in eth do 
		ll := Split(l,";");
		if ll[1] ne "level" then 
			p := StringToInteger(ll[2]);
			if p le 20 then 
				Append(~ethereal_ple20,l);
			end if;
			Append(~ethereal,l);
		end if;
	end for;
	return #ethereal, #ethereal_ple20;
end function;


for d in [1,2,3,7,11] do 
	eth, eth20 := EtherealData(d);
	lifts, lifts20 := AllLifts(d);
	ver, ver20 := VerifiedData(d);
	printf "$%o$ & $%o$ & $%o$ & $%o$ & $%o$ & $%o$ & $%o$ \\\\ \n", d, eth20, lifts20, ver20, eth, lifts, ver ;
end for;


