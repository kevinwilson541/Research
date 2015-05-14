#!/usr/bin/python

m = 2
p = 17
while not (p.is_prime() and (p-1) % (2*m)==0):
    p += 1
K = IntegerModRing(p)
q = (K.multiplicative_generator())^((p-1)/2/m)
print "(m,p,q)=",(m,p,q)

simple_reflection_names = "stu"
roots = ["alpha_"+ref for ref in simple_reflection_names]
R=PolynomialRing(K, len(roots), roots) # don't do this anymore

#I = R.base_ring().gens()[0]

C = Matrix([
    [   2,   -1,           -q**(2*m-1)],
    [  -1,    2,           -q],
    [  -q,   -q**(2*m-1),    2],
])

def alpha(x):
    return R.gens_dict()["alpha_"+str(x)]

simple_reflections=var(list(simple_reflection_names))

def a(i,j):
    assert i in simple_reflections, "%s isn't in simple reflections" % i
    assert j in simple_reflections, "%s isn't in simple reflections" % j
    u = simple_reflections.index(i)
    v = simple_reflections.index(j)
    return R(C[u,v])

EndR = Hom(R,R)

remember_reflection = {}
def reflection_in(i):
    if not (i in remember_reflection):
        remember_reflection[i] = EndR( [alpha(j) - a(i, j) * alpha(i) for j in simple_reflections] )
    return remember_reflection[i]

remember_demazure = {}
def demazure(x):
    if not (x in remember_demazure): 
        remember_demazure[x] = lambda f:R((f - reflection_in(x)(f))/alpha(x))
    return remember_demazure[x]


#================================== 
# Stuff above here should get reimplemented, with polynomials replaced by 
# big dictionaries of their values over R.

remember_fast_apply_list_demazure_monomial = {}
def fast_apply_list_demazure_monomial(L, m):
    key = (tuple(L), m)
    if key not in remember_fast_apply_list_demazure_monomial: 
        if len(L) == 1:
            result = demazure(L[0])(m)
        else:
            (prefix, suffix) = (L[:-1], L[-1])
            p = demazure(suffix)(m)
            result = apply_list_demazure(prefix, p)
        remember_fast_apply_list_demazure_monomial[key] = result 
    return remember_fast_apply_list_demazure_monomial[key]

def apply_list_demazure(L,p):
    if len(L) == 0:
        return p
    if p.parent() != R:
        p = R(p)
    coeffs = p.coefficients()
    monomials = p.monomials()
    result = R(0)
    for i in range(len(monomials)):
        result += coeffs[i] * fast_apply_list_demazure_monomial(L, monomials[i])
    return result
    
def slow_apply_list_demazure(L, p):
    result = p
    for x in reversed(L):
        result = demazure(x)(result)
    return result

# give me a list of all tuples in N^n that sum to k.
def points_in_simplex(n, k):
    if n == 1:
        return [(k,)]
    else:
        results = []
        for t in reversed(range(k+1)):
            for low_dim_point in points_in_simplex(n-1, k-t):
                p = list(low_dim_point)
                p.insert(0,t)
                results.append(tuple(p))
        return results

def monomials_of_degree(k):
    n = len(simple_reflections)
    results = []
    for multidegree in points_in_simplex(n,k):
        monomial = R(1)
        for i in range(n):
            var = simple_reflections[i]
            monomial *= alpha(var) ** multidegree[i]
        results.append(monomial)
    return results

def demazure_list(degree):
    if degree == 0:
        raise Exception("You don't ever need to think about degree 0")
    if degree == 1:
        return [[x] for x in list(simple_reflections)]
    shorter_list = demazure_list(degree-1)
    results = []
    for short_L in shorter_list:
        for x in simple_reflections:
            L = list(short_L)
            if x != L[-1]:
                L.append(x)
                results.append(L)
    return results

memoize_pruned_demazure_list = {}
def pruned_demazure_list(degree):
    if degree <= 2:
        return demazure_list(degree)   
    if degree not in memoize_pruned_demazure_list:
        shorter_list = pruned_demazure_list(degree-1)
        results = []
        for short_L in shorter_list:
            for x in simple_reflections:
                L = list(short_L)
                if x != L[-1]:
                    if x != L[-2]:
                        L.append(x)
                        results.append(L)
                    elif x==s:
                        L.append(x)
                        results.append(L)
                    elif x==t and L[-1]==u:
                        L.append(x)
                        results.append(L)
        memoize_pruned_demazure_list[degree] = results
    return memoize_pruned_demazure_list[degree] 

memoize_repruned_demazure_list = {}
def repruned_demazure_list(degree):
    if degree <= 4:
        return pruned_demazure_list(degree)   
    if degree not in memoize_repruned_demazure_list:
        shorter_list = repruned_demazure_list(degree-1)
        results = []
        for short_L in shorter_list:
            for x in simple_reflections:
                L = list(short_L)
                if x != L[-1]:
                    if x != L[-2]:
                        L.append(x)
                        results.append(L)
                    elif L[-1] == u and L[-4] == u:
                        pass
                    elif x==s:
                        L.append(x)
                        results.append(L)
                    elif x==t and L[-1]==u:
                        L.append(x)
                        results.append(L)
        memoize_repruned_demazure_list[degree] = results
    return memoize_repruned_demazure_list[degree] 


class NotReducedException(Exception):
    pass

suffix_repair = {
    (t,s,t):([s,t,s], 1), 
    (u,s,u):([s,u,s], q), 
    (u,t,u):([t,u,t], q^(2*m-1)), 
    (u,t,s,u,s):([t,u,t,s,u], q^(2*m-2)),
    (u,s,t,u,t):([s,u,s,t,u], q^2),
}

class SuppressedSuffixException(Exception):
    pass

def repair_suffix(word):
    for pattern in suffix_repair:    
        n = len(pattern)
        suffix = word[-n:]
        if tuple(suffix) == pattern:
            (repaired_suffix, const) = suffix_repair[pattern]
            return (word[:-n], repaired_suffix, const)
    raise SuppressedSuffixException("%s didn't have a suppressed suffix" % str(word))

memoize_append_element = {}
def append_element(reduced_expr, suffix_element):
    key = (tuple(reduced_expr), suffix_element)
    if key not in memoize_append_element:
        #print "I'm now appending %s to %s" % (suffix_element, reduced_expr)
        n = len(reduced_expr)
        longer_element = copy(reduced_expr)
        longer_element.append(suffix_element)
        #print "here", key
        if n == 0:
            memoize_append_element[key] = ([suffix_element], 1)
        else: 
            small_list = repruned_demazure_list(n)
            big_list = repruned_demazure_list(n+1)
            assert reduced_expr in small_list, "%s isn't in small list" % reduced_expr
            if suffix_element == reduced_expr[-1]:
                memoize_append_element[key] = (None, None)
                raise NotReducedException("%s is not reduced" % str(longer_element))
            if n == 1:
                assert longer_element in big_list, "%s isn't in big list" % str(longer_element)
                memoize_append_element[key]=(longer_element, 1)
            elif suffix_element == reduced_expr[-2]:
                try:
                    (prefix, preferred_suffix, suffix_const) = repair_suffix(longer_element)
                except SuppressedSuffixException:
                    assert longer_element in big_list, "%s isn't in big list" % str(longer_element)
                    memoize_append_element[key] = (longer_element, 1)
                else:
#                   This might well raise a NotReducedException, which we do not handle
                    try:
                        (new_word, prefix_const) = append_element(prefix, preferred_suffix[0])
                    except NotReducedException as E:
                        memoize_append_element[key] = (None, None)
                        raise E
                    else:
                        new_word.extend(preferred_suffix[1:])
                        assert new_word in big_list, "%s isn't in big list" % str(new_word)
                        memoize_append_element[key]=(new_word, prefix_const * suffix_const)
            else:
                assert longer_element in big_list, "%s isn't in big list" % str(longer_element)
                memoize_append_element[key]=(longer_element, 1)
    if memoize_append_element[key] == (None, None):
        raise NotReducedException()
    (result_word, result_const) = memoize_append_element[key] 
    return (copy(result_word), result_const)

def test_append_element(n):
    L = repruned_demazure_list(n)
    big_L = repruned_demazure_list(n+1)
    for word in L:
        for new_guy in [s,t,u]:
            try:
                (result, const) = append_element(word, new_guy)
                print word, "+", new_guy, "----->", result, const,
                candidate = copy(word)
                candidate.append(new_guy)
                if candidate != result:
                    print "**"
                else:
                    print
            except NotReducedException:
                print word, "+", new_guy, "is not reduced"


def normal_form_append(w, x):
    try:
        return append_element(w, x)
    except NotReducedException:
        return (None, None)

def test_normal_form(n):
    L = repruned_demazure_list(n)
    big_L = repruned_demazure_list(n+1)
    for word in L:
        for new_guy in [s,t,u]:
            (result, const) = normal_form_append(word, new_guy)
            print word, "+", new_guy, "----->", result, const,
            candidate = copy(word)
            candidate.append(new_guy)
            if candidate != result:
                print "**"
            else:
                print

def multiply_by(length, suffix):
    entries = {}
    L = repruned_demazure_list(length)        
    M = repruned_demazure_list(length+1)
    M_lookup={}
    for j in range(len(M)):
        M_lookup[tuple(M[j])] = j

    for i in range(len(L)):
        demazure = L[i]
        (result, coeff) = normal_form_append(demazure, suffix)
        if result != None:
            j = M_lookup[tuple(result)]
            entries[(i, j)] = coeff
    return Matrix(len(L), len(M), entries)


def find_all_relation_vectors(degree_of_demazure, degree_of_input_monomial):
    results = []
    V = None
    for monomial in monomials_of_degree(degree_of_input_monomial):
        row = []
        for operators in repruned_demazure_list(degree_of_demazure):
            result = fast_apply_list_demazure_monomial(operators, monomial)
            row.append(result)
        results.append(row)
    for output_monomial in monomials_of_degree(degree_of_input_monomial-degree_of_demazure):
        coeff_matrix = []
        for raw_row in results:
            row = []
            for poly in raw_row:
                row.append(K(poly.coefficient(output_monomial)))
            coeff_matrix.append(row)
        new_nullspace = matrix(coeff_matrix).right_kernel()
        #print new_nullspace
        if V==None:
            V = new_nullspace
        else:
            V = V.intersection(new_nullspace)
    return V

#Assume that V is the output of find_all_relation_vectors(old_degree_of_demazure, BIG)
def find_boring_relation_vectors(V, old_degree_of_demazure):
    S = linear_transformation(multiply_by(old_degree_of_demazure, s))
    T = linear_transformation(multiply_by(old_degree_of_demazure, t))
    U = linear_transformation(multiply_by(old_degree_of_demazure, u))
    return span(S(V).basis() + T(V).basis() + U(V).basis())

A = FreeAlgebra(K, 3, simple_reflections)
d = {
    s:A.gens()[0],   t:A.gens()[1],   u:A.gens()[2],
    's':A.gens()[0], 't':A.gens()[1], 'u':A.gens()[2]
}
f = lambda L:prod([d[x] for x in L])

def dot_product(v, L):
   return sum([x*y for (x,y) in zip(v, L)])

def show_relations(V, demazure_length, min_relation_length=0):
    L = [f(x) for x in repruned_demazure_list(demazure_length)]
    results = []
    for relation_vec in V.basis():
        relation = dot_product(relation_vec, L)
        if len(relation.terms()) > min_relation_length:
            results.append(relation) 
    return results


#def find_all_relations(degree_of_demazure, degree_of_input_monomial, min_relation_length=0):
#    V = find_all_relation_vectors(degree_of_demazure, degree_of_input_monomial)
#    return show_relations(V, degree_of_demazure, min_relation_length)


# both parameters vector spaces, please
def find_interesting_relations(boring_relations, all_relations):
    boring_count = boring_relations.dimension()
    M_boring = boring_relations.basis_matrix()
    M_all = all_relations.basis_matrix()
    M = M_boring.stack(M_all).transpose()
    pivot_list = M.pivots()
    for i in range(boring_count):
        assert pivot_list[i] == i 
    interesting_pivots = pivot_list[boring_count:]
    interesting_matrix = M.matrix_from_columns(interesting_pivots).transpose()
    return all_relations.subspace(interesting_matrix.rows())


def find_all_relations_stable(degree_of_demazure):
    monomial_degree = degree_of_demazure - 1
    V = None
    happiness = 0
    while(happiness < 2): 
        monomial_degree += 1
        print "acting on monomials of degree", monomial_degree
        relation_space = find_all_relation_vectors(degree_of_demazure, monomial_degree)
        if V == None:
            V = relation_space
        else:
            if V.is_subspace(relation_space):
                happiness += 1
            else:
                V = V.intersection(relation_space)
                happiness = 0 #damnit
    return V

def classify_relations():
    degree_of_demazure = 0
    results = {}
    relations = VectorSpace(QQ, 0)     
    poincare_coeffs = [1] # in degree 0, we have the identity

    while relations.dimension() == 0:
        degree_of_demazure += 1
        print "Trying Demazure degree", degree_of_demazure, "..."
        relations = find_all_relations_stable(degree_of_demazure)
        results[degree_of_demazure] = relations
        print "Found %d relations." % relations.dimension()
        coeff = len(repruned_demazure_list(degree_of_demazure))-relations.dimension()
        poincare_coeffs.append(coeff)
        print "Poincare coefficient", coeff

    print "Okay, now we're getting somewhere. Relations are:"
    print show_relations(relations, degree_of_demazure)
    new_relation_count = relations.dimension()
    while (poincare_coeffs[-1] > 0):
        degree_of_demazure += 1
        print "Trying Demazure degree", degree_of_demazure, "..."
        old_relations = relations
        relations = find_all_relations_stable(degree_of_demazure)
        boring = find_boring_relation_vectors(old_relations, degree_of_demazure-1)
        interesting = find_interesting_relations(boring, relations)
        results[degree_of_demazure] = interesting
        new_relation_count = interesting.dimension()
        print "Found %d interesting relations, as follows:" % interesting.dimension()
        print show_relations(interesting, degree_of_demazure)
        coeff = len(repruned_demazure_list(degree_of_demazure))-relations.dimension()
        poincare_coeffs.append(coeff)
        print "Poincare coefficient", coeff
         
    return (results, poincare_coeffs)

def check_relation(relation_vec, degree_of_demazure, degree_of_monomials):
    L = repruned_demazure_list(degree_of_demazure)
    for monomial in monomials_of_degree(degree_of_monomials):
        result = R(0)
        for j in range(len(L)):
            result += relation_vec[j] * apply_list_demazure(L[j], monomial)
        print result
