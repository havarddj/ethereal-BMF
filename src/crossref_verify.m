
AttachSpec("BMF/spec");

// cross-referencing verified data with the nonEis data files 
procedure CrossReferenceVerify(d)

	K := QuadFld(d);
	ZK := MaximalOrder(K);

	// load in the verified and nonEis data 
	printf "Loading verified data for d = %o...\n", d;
	verified := Split(Read("data/verified_d" cat Sprint(d) cat ".csv"));

	printf "Loading non-Eisenstein data for d = %o...\n", d;
	nonEis := Split(Read("data/nonEis_d" cat Sprint(d) cat ".csv"));

	// these are the actual ideals represented by our LMFDB label strings 
	nn := #Split(nonEis[1],";");
	idealStrings := Split(nonEis[1],";")[3..nn-3];
	ideals := [LMFDBIdeal(K,u) : u in idealStrings];


	outFileName := "data/new_verified_d" cat Sprint(d) cat ".csv";

	// checks if the verified data ver matches the data non from nonEis
	// only checks prime not over the level and characteristic! 
	function IsMatch(ver,eigSeq,non)
		sver := Split(ver,";");
		snon := Split(non,";");

		level := LMFDBIdeal(K,sver[1]);
		p := StringToInteger(sver[2]);

		match := true;
		for i in [1..#ideals] do 
			if GCD(ideals[i],level*p) eq 1*ZK then 
				match and:= snon[2+i] eq Sprint(eigSeq[i]);
			else
				//print LMFDBLabel(ideals[i]);
			end if;
		end for;

		return match;
	end function;

	// conjugate of ideal 
	function ConjugateIdeal(id)
		a := Automorphisms(K)[2];
		return ideal<ZK | [a(u) : u in Basis(id)]>;
	end function;

	// conjugates the label of the ideal 
	function ConjugateIdealLabel(lab)
		return LMFDBLabel(ConjugateIdeal(LMFDBIdeal(K,lab)));
	end function;


	// performs the complex conjugation operation on the whole verified 
	// string. this changes the level to its conjugate, and swaps eigenvalues 
	// over split primes
	function ConjugateVerData(ver,eigSeq)
		sver := Split(ver,";");
		conjLevel := ConjugateIdealLabel(sver[1]);
		conjLiftLabel := ConjugateIdealLabel(sver[3]);
		conjIdealStrings := [ConjugateIdealLabel(u) : u in idealStrings];

		// sadly can't pass parallel sort our own ordering function
		// and the default one doesn't sort labels how we want! 
		conjEigStr := [eigSeq[Index(conjIdealStrings,u)] : u in idealStrings];

		// and to match the rest of the verified data, we need Sprint(conjEigStr) without spaces!
		s := Sprint(conjEigStr);
		ss := "[" cat s[3..#s-2] cat "]";

		conjver := conjLevel cat ";" cat sver[2] cat ";" cat conjLiftLabel cat ";" cat ss;
		return conjver,conjEigStr;
	end function;

	// now we go through the verified data one row at a time and try to match it to 
	// something in nonEis 
	conjCount := 0;
	SetColumns(0);
	print "Checking verified data against non-Eisenstein...";
	problems := false;
	for i in [1..#verified] do 
		v := verified[i];
		sv := Split(v,";");
		level := LMFDBIdeal(K,sv[1]);
		eigSeq := [StringToInteger(u) : u in Split(sv[4][2..#sv[4]-1],",")];


		// we collect everything with the same level norm and prime. this should be a short list! 
		possible := [];
		for c in nonEis[2..#nonEis] do 
			ee := Split(c,";");
			if Norm(LMFDBIdeal(K,ee[1])) eq Norm(level) and ee[2] eq sv[2] then 
				Append(~possible,c);
			end if;
		end for;

		// the conjugated data and the entries in `possible` which match (should be one per)
		conjv,conjEigSeq := ConjugateVerData(v,eigSeq);
		match_v := [IsMatch(v,eigSeq,u) : u in possible];
		match_conjv := [IsMatch(conjv,conjEigSeq,u) : u in possible];

		// we get hold of the matching system 
		if true in match_v then
			match := possible[Index(match_v,true)];
			ss := Split(v,";");
		elif true in match_conjv then 
			match := possible[Index(match_conjv,true)];
			ss := Split(conjv,";");
			conjCount +:= 1;
		else 
			printf "PROBLEM: the following verified entry does not apper in nonEis: %o\n", v;
			problems := true;
			if Split(v,";")[1] eq Split(v,";")[3] then 
				print "Lift appears to be trivial, i.e. form was never ethereal";
			end if;
			print "";
		end if;


		// we get the eigenvalues from the match, and write these instead (since they contain the 
		// bad prime data that was previously missing)
		mm := Split(match,";");
		eigs := Sprint(mm[3..#mm-3]);
		out_v := ss[1] cat ";" cat ss[2] cat ";" cat ss[3] cat ";" cat "[" cat eigs[3..#eigs-2] cat "]";
		Write(outFileName,out_v);
	end for;

	printf "Conjugated %o forms to minimal LMFDB ideal.\n",conjCount;

	if problems then
		print "Problems occured during this run. Some verified systems may not have been matched to existing data.";
	else 
		print "No problems occured";
	end if;
	print "";

end procedure;


for d in [1,2,3,7,11] do 
	CrossReferenceVerify(d);
end for;

