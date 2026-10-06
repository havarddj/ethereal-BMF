

AttachSpec("../BMF/spec");

// helper function for turning magma printouts into LaTeX strings
// str1 is the string Magma uses for an element in a field
// str2 is the string we want to replace it with 
function RenderEig(eig,str1,str2)
	eigstr := Sprint(eig);

	if str1 ne str2 then 
		// first we replace occurences of str1 with str2 
		while str1 in eigstr do 
			ind := Index(eigstr,str1);
			eigstr := eigstr[1..ind-1] cat str2 cat eigstr[ind+#str1..#eigstr];
		end while;
	end if;

	// then we delete all asterisks
	while "*" in eigstr do
		ind := Index(eigstr,"*");
		eigstr := eigstr[1..ind-1] cat eigstr[ind+1..#eigstr];
	end while;

	return eigstr;
end function;


// turns the sequence `seq` into a LaTeX table row 
function TableRow(seq)
	row := "";
	for i in [1..#seq] do 
		row cat:= seq[i];
		if i ne #seq then 
			row cat:= " & ";
		else 
			row cat:= " \\\\";
		end if;
	end for;
	return row;
end function;



procedure Table_1()

	d := 1;
	K := QuadFld(d);
	ZK := MaximalOrder(K);

	F<w> := NumberField(Polynomial([-1,1,1]));
	ZF := MaximalOrder(F);
	FF<wb>, qq := ResidueClassField(3*ZF);
	SetPowerPrinting(FF,false);

	level := LMFDBIdeal(K,"293.1");

	print "Computing mod p cohomology";
	B := BianchiCohomologySpace(level, BianchiWeight(K,0,0 : char := 3));
	BB := ChangeRing(B,FF);

	labels := ["2.1", "5.1", "9.1", "13.1", "13.2", "17.1", "17.2", "29.1", "29.2", "37.1"];
	PP := [LMFDBIdeal(K,l) : l in labels];

	// pick out the ethereal form 
	HP := [HeckeOperator(BB,P) : P in PP];
	Ep := Eigenforms(BB);
	EEp := [e : e in Ep | Eigenvalue(e,PP[1]) eq (2*wb + 2)];
	assert #EEp eq 1;
	f := EEp[1];

	print "Computing char. 0 cohomology";
	// and the characteristic 0 lift 
	B0 := BianchiCohomologySpace(level * LMFDBIdeal(K,"5.2"), BianchiWeight(K,0,0));
	H0 := [HeckeOperator(B0,P) : P in PP];
	E0 := Eigenforms(B0);

	_,mm := IsIsomorphic(CoefficientRing(E0[1]),F);

	EE0 := [e : e in E0 | mm(Eigenvalue(e,PP[1])) eq -w - 1];
	assert #EE0 eq 1;
	f0 := EE0[1];

	
	SetColumns(0);
	// we build all four rows simultaneously 
	toprow := "$\\p$ & ";
	row1 := "$a_{\\p}(f)$ & ";
	row2 := "$a_{\\p}(F) \\Mod{3}$ & ";
	row3 := "$a_{\\p}(F)$ & ";
	for l in labels do 
		toprow cat:= "\\texttt{" cat l cat "}";
		P := LMFDBIdeal(K,l);
		row1 cat:= "$" cat RenderEig(Eigenvalue(f,P),"wb","\\bar{w}") cat "$";
		row2 cat:= "$" cat RenderEig(qq(mm(Eigenvalue(f0,P))),"wb","\\bar{w}") cat "$";
		row3 cat:= "$" cat RenderEig(mm(Eigenvalue(f0,P)),"","") cat "$";

		if l ne labels[#labels] then 
			toprow cat:= " & ";
			row1 cat:= " & ";
			row2 cat:= " & ";
			row3 cat:= " & ";
		else 
			toprow cat:= "\\\\ \n\\midrule";
			row1 cat:= "\\\\";
			row2 cat:= "\\\\";
			row3 cat:= "\\\\";
		end if;
	end for;

	print "Printing contents of table 1:";
	print "";
	print toprow;
	print row1;
	print row2;
	print row3;

end procedure;




// table 2

d := 3;
K := QuadFld(d);
ZK := MaximalOrder(K);

F<w> := NumberField(Polynomial([-2,0,1]));
ZF := MaximalOrder(F);
P := ideal<ZF | [1,4]>;
FF<wb>, qq := ResidueClassField(P);

level := LMFDBIdeal(K,"607.1");

print "Computing mod p cohomology";
B := BianchiCohomologySpace(level, BianchiWeight(K,0,0 : char := 31));

labels := ["3.1", "4.1", "7.1", "7.2", "13.1", "19.1", "19.2", "25.1", "31.1", "31.2"];
PP := [LMFDBIdeal(K,l) : l in labels];

// pick out the ethereal form 
HP := [HeckeOperator(B,P) : P in PP];
Ep := Eigenforms(B);
EEp := [e : e in Ep | Eigenvalue(e,PP[1]) eq 16];
assert #EEp eq 1;
f := EEp[1];

print "Computing char. 0 cohomology";
// and the characteristic 0 lift 
B0 := BianchiCohomologySpace(level * LMFDBIdeal(K,"13.2"), BianchiWeight(K,0,0));
H0 := [HeckeOperator(B0,P) : P in PP];
E0 := Eigenforms(B0);

_,mm := IsIsomorphic(CoefficientRing(E0[1]),F);

EE0 := [e : e in E0 | mm(Eigenvalue(e,PP[1])) eq -w - 1];
assert #EE0 eq 1;
f0 := EE0[1];





procedure Table_3()
	print "Printing contents of table 3:";
	print "$d$& $2, 7, 11$ & $1$ & $3$  \\\\ \\midrule";
    print "$\\norm(\\mf{n})$& $50{,}000$ & $100{,}000$ & $150{,}000$ \\\\";
end procedure;



// table 4 
d := 7;
K := QuadFld(d);
ZK := MaximalOrder(K);

F<w> := NumberField(Polynomial([1,-9,8,17,-6,-8,1,1]));
ZF := MaximalOrder(F);
P := ideal<ZF | w^6 - 8*w^4 + w^3 + 14*w^2 - 2>;
assert Norm(P) eq 53;

FF,qq := ResidueClassField(P);

level := LMFDBIdeal(K,"667.2");
B := BianchiCohomologySpace(level, BianchiWeight(K,0,0 : char := 53));

labels := ["2.1", "2.2", "7.1", "9.1", "11.1", "11.2", "23.1", "23.2", "25.1", "53.1", "53.2"];
PP := [LMFDBIdeal(K,l) : l in labels];

liftLevel := level * LMFDBIdeal(K,"7.1");
B0 := BianchiCohomologySpace(liftLevel, BianchiWeight(K,0,0));
H0 := [HeckeOperator(B0,P) : P in PP];
E0 := Eigenforms(B0);

EE0 := [u : u in E0 | Eigenvalue(u,PP[1]) eq 4 - 10*w - 2*w^2 + 7*w^3 - w^5];
assert #EE0 eq 1;
f0 := EE0[1];


toprow := "$\\p$ & $a_{\\p}(f)$ & $a_{\\p}(F) \\Mod{ \\mf{p}_{53}}$ & $a_{\\p}(F)$ \\\\ \\midrule";
for l in labels do 
	row := l;
	P := LMFDBIdeal(K,l);
	if GCD(P,liftLevel*53) ne 1*ZK then 
		row cat:= "${}^{\\dagger}$";
	end if;

	row cat:= " & " cat Sprint(Eigenvalue(f,P));

	if GCD(P,liftLevel) ne 1*ZK then 
		row cat:= "$-$ & $-$";
	else 
		row cat:= "$" cat Sprint(qq(Eigenvalue(f0,P))) cat "$ & $" cat RenderEig(Eigenvalue(f0,P),"","") cat "$";
	end if;
	print row;

end for;




// table 5 

// table 6

// table 7

procedure Table_7()
	d := 7;
	K<w> := QuadFld(d);
	ZK := MaximalOrder(K);
	level := LMFDBIdeal(K,"512.1");

	B := BianchiCohomologySpace(level, BianchiWeight(K,0,0 : char := 7));

	labels := ["2.2", "9.1", "11.1", "11.2", "23.1", "23.2", "25.1", "29.1", "29.2", "37.1", "37.2"];
	PP := [LMFDBIdeal(K,l) : l in labels];
	HP := [HeckeOperator(B,P) : P in PP];
	E := Eigenforms(B);
	EE := [e : e in E | Eigenvalue(e,PP[3]) eq 6];
	assert #EE eq 1;
	f := EE[1];

	L := ext<K | Polynomial([w,0,1])>;
	ZL := MaximalOrder(L);
	decomps := [DecompositionType(ZL,P) : P in PP];

	tr := ["$\\p$"] cat ["\\texttt{" cat l cat "}" : l in labels];
	r1 := ["$\\mf{p}\\OO_L$"] cat [#dd eq 2 select "split" else "inert" : dd in decomps];
	r2 := ["$a_{\\p}(f)$"] cat ["$" cat Sprint(Eigenvalue(f,P)) cat "$" : P in PP];

	toprow := TableRow(tr) cat "\\midrule";
	row1 := TableRow(r1);
	row2 := TableRow(r2);


	print "Printing contents of table 7:";
	print toprow;
	print row1;
	print row2;
end procedure;


procedure Table_9()
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

	print "Printing contents of table 9:";
	print "$d$  & \\shortstack{Ethereal \\\\ ($p < 20$)} & Lifts & Verified &  \\shortstack{Ethereal \\\\ (all $p$)}& Lifts & Verified  \\\\ \\midrule";
	for d in [1,2,3,7,11] do 
		eth, eth20 := EtherealData(d);
		lifts, lifts20 := AllLifts(d);
		ver, ver20 := VerifiedData(d);
		assert &and([checklifts(u) : u in [lifts,lifts20]] cat [checkver(u) : u in [ver,ver20]]);
		printf "$%o$ & $%o$ & $%o$ & $%o$ & $%o$ & $%o$ & $%o$ \\\\ \n", d, #eth20, #lifts20, #ver20, #eth, #lifts, #ver ;
	end for;
	print "";

end procedure;