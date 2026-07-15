
load "look_for_lifts.m";

r := Read("data/nonEis_d11.csv");

lines := Split(r,"\n");

problems := [];


N := #Split(lines[1],";");

for l in lines do 
	ll := Split(l,";");

	if #ll ne N then 
		print "problem!!!";
		l;
		Append(~problems,<ll[1],StringToInteger(ll[2])>);
	end if;
end for;




