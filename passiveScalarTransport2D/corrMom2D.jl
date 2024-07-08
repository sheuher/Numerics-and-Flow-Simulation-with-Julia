function corrMom2D!(rhou, rhouP, p, dt, delta)
	# corrects momentum prediction field based on given pressure field
	# use CDS for pressure gradient

	nG = 1
	ImaAll, JmaAll = size(rhou)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	rhou[Ifi:Ila, Jfi:Jla, 1] = rhouP[Ifi:Ila, Jfi:Jla, 1] .- dt * 0.5 /delta *(p[Ifip:Ilap, Jfi:Jla] .-p[Ifim:Ilam, Jfi:Jla]) 

	rhou[Ifi:Ila, Jfi:Jla, 2] = rhouP[Ifi:Ila, Jfi:Jla, 2] .- dt * 0.5 /delta *(p[Ifi:Ila, Jfip:Jlap] .-p[Ifi:Ila, Jfim:Jlam]) 
	
	return nothing
end