from esutil.coords import randsphere

ra_min = 0. 
dec_min = 0.
ra_max = 90.
dec_max = 90.
n_rand = 1000.
ra,dec = randsphere(n_rand, ra_range=[ra_min,ra_max], dec_range=[dec_min, dec_max])
print ra, dec

