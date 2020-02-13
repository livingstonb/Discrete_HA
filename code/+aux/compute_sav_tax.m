function sav_tax = compute_sav_tax(sav, taxrate, thresh)
	taxable_sav = max(sav - thresh, 0);
	sav_tax = taxrate * taxable_sav;
end