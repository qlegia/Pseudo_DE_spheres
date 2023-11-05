load FHY100s;
soln_noh = sum(FHYnoh);
goodsoln_noh = sum(g_FHYnoh);
err_noh = max(abs(soln_noh-goodsoln_noh));
for q=1:max_q
   soln(q,:) = sum(FHY{q}); 
      % sum up all ell's and m's at each point for the solution
   goodsol(q,:) = sum(goodFHY{q});
   err(q) = max(abs(soln(q,:)-goodsol(q,:)));
   exact_err(q) = max(abs(goodsol(q,:)));
   exact_err2(q) = max(abs(soln(q,:)));
end  

