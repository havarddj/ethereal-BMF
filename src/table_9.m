
AttachSpec("../BMF/spec");


function VerifiedData(d)
	K := QuadFld(d);
	ZK := MaximalOrder(K);

	input := "../data/verified_d" cat Sprint(d) cat ".csv";
	lines := Split(Read(input), "\n");

	p_le20 := [];
	for l in lines do 
		ll := Split(l,";");
		p := StringToInteger(ll[2]);
		if p le 20 then 
			Append(~p_le20,l);
		end if;
	end for;
	return lines, p_le20;
end function;


function AllLifts(d)
	K := QuadFld(d);
	ZK := MaximalOrder(K);

	lmfdb_lifts := Split(Read("../sage/d" cat Sprint(d) cat "_liftable.csv"));
	lmfdb_lifts := lmfdb_lifts[2..#lmfdb_lifts];
	irrat_data := Split(Read("../data/irrat_lifts_d" cat Sprint(d) cat ".csv"));

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
			if ll[3] ne "None found" and ll[3] ne "No level raising primes available" then 
				p := StringToInteger(ll[2]);
				if p le 20 then 
					Append(~irrat_lifts_ple20,l);
				end if;
				Append(~irrat_lifts,l);
			end if;
		end if;
	end for;

	return lmfdb_lifts cat irrat_lifts, lmfdb_lifts_ple20 cat irrat_lifts_ple20;
end function;


function EtherealData(d)
	K := QuadFld(d);
	ZK := MaximalOrder(K);

	eth := Split(Read("../data/nonEis_d" cat Sprint(d) cat ".csv"));

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
	return ethereal, ethereal_ple20;
end function;


foo := function(str)
	ss := Split(str,";");
	return ss[1] cat ";" cat ss[2] cat ";" cat ss[3];
end function;

bar := function(str)
	ss := Split(str,";");
	return ss[1] cat ";" cat ss[2] cat ";" cat ss[4];
end function;

checklifts := function(list)
	ll := list;
	return &and [#[u : u in ll | u eq ll[i]] eq 1 : i in [1..#ll]];
end function;

checkver := function(list)
	ll := [bar(u) : u in list];
	return &and [#[u : u in ll | u eq ll[i]] eq 1 : i in [1..#ll]];
end function;


for d in [1,2,3,7,11] do 
	eth, eth20 := EtherealData(d);
	lifts, lifts20 := AllLifts(d);
	ver, ver20 := VerifiedData(d);
	assert &and([checklifts(u) : u in [lifts,lifts20]] cat [checkver(u) : u in [ver,ver20]]);
	printf "$%o$ & $%o$ & $%o$ & $%o$ & $%o$ & $%o$ & $%o$ \\\\ \n", d, #eth20, #lifts20, #ver20, #eth, #lifts, #ver ;
end for;



d := 2;
eth, eth20 := EtherealData(d);
lifts, lifts20 := AllLifts(d);
ver, ver20 := VerifiedData(d);

