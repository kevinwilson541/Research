#!/usr/bin/python

class ModPoly:
    def __init__(self, ring, coeffs=dict([((0,0,0), 0)]), poly=None, variables=var('x,y,z')):
        self.__x, self.__y, self.__z = variables        
        self.__ring = ring.base_ring()
        self.__polyRing = ring
        self.__valTable = dict()
        if poly:
            poly = self.__polyRing(poly)
            self.from_polynomial(poly)
        else:
            self.__coeffs = coeffs
            self.__generateValues(self.__valTable, self.__coeffs)
            self.partials = {
                x: dict(),
                y: dict(),
                z: dict()
            }
            
            for elem in self.partials:
                self.__generatePartial(elem)

    def base_ring(self):
        return self.__ring

    def parent(self):
        return self.__polyRing

    def coefficients(self):
        return self.__coeffs

    def table(self):
        return self.__valTable

    def variables(self):
        return (self.__x, self.__y, self.__z)
    
    def __generateValues(self, valTable, coeffs):
        for i in self.__ring:
            for j in self.__ring:
                for k in self.__ring:
                    total = 0
                    for elem in coeffs:
                        total += coeffs[elem]*i**elem[0]*j**elem[1]*k**elem[2]
                    valTable[(i,j,k)] = total

    def __generatePartial(self, var):
        coeffs = dict()
        for elem in self.__coeffs:
            if var == x:
                if elem[0] != 0: 
                    coeffs[(elem[0]-1, elem[1], elem[2])] = elem[0]*self.__coeffs[elem]
            elif var == y:
                if elem[1] != 0:
                    coeffs[(elem[0],elem[1]-1, elem[2])] = elem[1]*self.__coeffs[elem]
            else:
                if elem[2] != 0:
                    coeffs[(elem[0], elem[1], elem[2]-1)] = elem[2]*self.__coeffs[elem]

        self.__generateValues(self.partials[var], coeffs)

    def is_zero(self):
        return reduce(lambda x, y: x and self.__valTable[y].is_zero(), self.__valTable, True)

    def __add__(self,rhs):
        assert self.__polyRing == rhs.__polyRing
        coeffs = dict()
        for elem in self.__coeffs:
            if rhs.__coeffs.has_key(elem):
                coeffs[elem] = self.__coeffs[elem]+rhs.__coeffs[elem]
            else:
                coeffs[elem] = self.__coeffs[elem]
        for elem in rhs.__coeffs:
            if not coeffs.has_key(elem):
                coeffs[elem] = rhs.__coeffs[elem]
        return ModPoly(self.__polyRing, coeffs)

    def __sub__(self, rhs):
        assert self.__polyRing == rhs.__polyRing
        coeffs = dict()
        for elem in self.__coeffs:
            if rhs.__coeffs.has_key(elem):
                coeffs[elem] = self.__coeffs[elem] - rhs.__coeffs[elem]
            else:
                coeffs[elem] = self.__coeffs[elem]
        for elem in rhs.__coeffs:
            if not coeffs.has_key(elem):
                coeffs[elem] = -1*rhs.__coeffs[elem]
        return ModPoly(self.__polyRing, coeffs)

    def __mul__(self, rhs):
        assert rhs in self.__ring
        coeffs = dict()
        for elem in self.__coeffs:
            coeffs[elem] = self.__coeffs[elem]*rhs
        return ModPoly(self.__polyRing, coeffs)

    def __div__(self, rhs):
        assert rhs in self.__ring
        coeffs = dict()
        for elem in self.__coeffs:
            coeffs[elem] = self.__coeffs[elem]/rhs
        return ModPoly(self.__polyRing, coeffs)

    def __eq__(self, rhs):
        assert self.__polyRing == rhs.__polyRing
        return self.__valTable == rhs.__valTable 

    def __ne__(self, rhs):
        assert self.__polyRing == rhs.__polyRing
        return self.__valTable != rhs.__valTable

    def __call__(self, x, y, z):
        assert x in self.__ring
        assert y in self.__ring
        assert z in self.__ring
        return self.__valTable[(x,y,z)] 

    def __iadd__(self, rhs):
        assert self.__polyRing == rhs.__polyRing
        for elem in self.__coeffs:
            if rhs.__coeffs.has_key(elem):
                self.__coeffs[elem] += rhs.__coeffs[elem]
        for elem in rhs.__coeffs:
            if not self.__coeffs.has_key(elem):
                self.__coeffs[elem] = rhs.__coeffs[elem]
        
        self.__valTable = dict()
        self.__generateValues(self.__valTable, self.__coeffs)
        for elem in self.partials:
            self.__generatePartial(elem)
        return self

    def __isub__(self, rhs):
        assert self.__polyRing == rhs.__polyRing
        for elem in self.__coeffs:
            if rhs.__coeffs.has_key(elem):
                self.__coeffs[elem] -= rhs.__coeffs[elem]
        for elem in rhs.__coeffs:
            if not self.__coeffs.has_key(elem):
                self.__coeffs[elem] = -1*rhs.__coeffs[elem]
        
        self.__valTable = dict()
        self.__generateValues(self.__valTable, self.__coeffs)
        for elem in self.partials:
            self.__generatePartial(elem)
        return self

    def __imul__(self, rhs):
        assert rhs in self.__ring
        for elem in self.__coeffs:
            self.__coeffs[elem] *= rhs

        self.__valTable = dict()
        self.__generateValues(self.__valTable, self.__coeffs)
        for elem in self.partials:
            self.__generatePartial(elem)
        return self

    def __idiv__(self, rhs):
        assert rhs in self.__ring
        for elem in self.__coeffs:
            self.__coeffs[elem] /= rhs

        self.__valTable = dict()
        self.__generateValues(self.__valTable, self.__coeffs)
        for elem in self.partials:
            self.__generatePartial(elem)
        return self

    def polynomial(self):
        return self.__polyRing(reduce(lambda prev, curr: prev + self.__coeffs[curr]*self.__x**curr[0]*self.__y**curr[1]*self.__z**curr[2], self.__coeffs, 0))

    def from_polynomial(self, poly):
        assert len(poly.parent().variable_names()) == 3
        assert getattr(poly, 'dict', None) != None
        self.__polyRing = poly.parent()
        self.__ring = poly.base_ring()
        poly = self.__polyRing(poly)

        self.__x, self.__y, self.__z = var(self.__polyRing.variable_names())

        self.__coeffs = poly.dict()
        self.__generateValues(self.__valTable, self.__coeffs)

        self.partials = {
            x: dict(),
            y: dict(),
            z: dict()
        }
        for elem in self.partials:
            self.__generatePartial(elem)
        
    def divX(self):
        print "stuff"

    def divY(self):
        print "stuff"

    def divZ(self):
        print "stuff"

