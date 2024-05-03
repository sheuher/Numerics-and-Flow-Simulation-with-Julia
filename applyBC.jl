function applyBC!(Phi)
    """
        apply periodic boundary condition
    """
    nG = 1
    ImaAll, JmaAll = size(Phi)
	Ifi = nG+1; Ifip = Ifi+1; Ifim = Ifi-1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfip = Jfi+1; Jfim = Jfi-1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1
    
    Phi[Ifim, Jfi:Jla] .= Phi[Ila, Jfi:Jla]; Phi[Ilap, Jfi:Jla] .= Phi[Ifi, Jfi:Jla]
    Phi[Ifi:Ila, Jfim] .= Phi[Ifi:Ila, Jla]; Phi[Ifi:Ila, Jlap] .= Phi[Ifi:Ila, Jfi]
end