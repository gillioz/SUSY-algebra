# SUSY algebra module

# use the SymPy library for manipulation of symbolic objects
from sympy import Rational  # rational numbers
from sympy import simplify  # simplification of symbolic expressions
from sympy import symbols   # symbolic expressions
from sympy import linsolve  # solving linear systems of equations
from sympy import sqrt      # square root
from sympy import Matrix    # matrix
from sympy import diff      # derivative
from sympy import I, re, im # complex numbers
from sympy import init_printing # for LaTeX output
from sympy.combinatorics import Permutation # permutations
import datetime             # for functions displaying results

# if the module IPython is available, then use it for display in latex:
#import importlib
#IPython_found = importlib.util.find_spec("IPython")
#if IPython_found is not None:
#    from IPython import display
#else:
#    def display(expression):
#        return print(expression)

init_printing()

# fractions defined symbolically
half = Rational(1,2)
quarter = Rational(1,4)

# global variable defining the type of output
latexoutput = False

# Minkowski metric
eta = ((1, 0, 0, 0),
       (0, -1, 0, 0),
       (0, 0, -1, 0),
       (0, 0, 0, -1))
# Pauli matrices
sigma = ( ((1, 0), (0, 1)),
          ((0, 1), (1, 0)),
          ((0, -I), (I, 0)),
          ((1, 0), (0, -1)) )
sigmabar = ( ((1, 0), (0, 1)),
             ((0, -1), (-1, 0)),
             ((0, I), (-I, 0)),
             ((-1, 0), (0, 1)) )
# delta and epsilon tensors
delta = ( (1, 0), (0, 1) )
epsilon = ( (0, 1), (-1, 0) )
# Levi-Civita tensor
def LeviCivita(indices):
    if len(set(indices)) < len(indices):
        return 0
    return Permutation(indices).signature()

# Generates a table of n indices
# taking values in the range [0, r)
# (the table will containt r^n elements)
def sumindices(r, n):
    return [ [ (i // r**j) % r for j in range(n)]
             for i in range(r**n) ]


# function used to evaluate a sum:
# for all values of the spin and Lorentz indices,
# prefactor is evaluated first, and then
# expression is evaluated only if prefactor is non-zero
def sum(prefactor, expression, nundotted, ndotted, nlorentz):
    """Evaluates the sum of prefactor * expression over a number of indices
    for all values of the spin and Lorentz indices,
    prefactor is evaluated first, and then
    expression is evaluated only once if prefactor is non-zero
    """
    # all possible summation indices
    undottedindices = sumindices(2, nundotted)
    dottedindices = sumindices(2, ndotted)
    lorentzindices = sumindices(4, nlorentz)
    first = True
    result = 0
    # find the cases for which the prefactor is non-zero
    for a in undottedindices:
        for b in dottedindices:
            for mu in lorentzindices:
                x = prefactor(a, b, mu)
                if x != 0: # evaluate expression and add the result
                    if first:
                        result = x * expression(a, b, mu)
                        first = False
                    else:
                        result += x * expression(a, b, mu)
    return result


##############################################################################
class Error(Exception):
    """Exceptions
    """
    def __init__(self, fct, message, expression):
        self.fct = fct
        self.message = message
        self.expression = expression


##############################################################################

class SingleField:
    """Object with a space-time dependence
    """

    def __init__(self, name, s1 = [], s2 = [], dim = 0):
        """Constructor for the class SingleField
        """
        self.name = name    # name
        self.s1 = s1[:]     # undotted spin indices
        self.s2 = s2[:]     # dotted spin indices
        self.der = []       # derivatives acting on field
        self.Grassmann = (len(s1) + len(s2)) % 2 == 1
                            # Grassmann number or not
        self.dim = dim      # scaling dimension

    def copy(self):
        """Creates a copy of the field
        """
        other = SingleField(self.name, self.s1, self.s2, self.dim)
        other.der = self.der[:]
        other.Grassmann = self.Grassmann
        return other

    def d(self, mu):
        """Defines how the derivative acts on a field
        """
        other = self.copy()
        other.der.append(mu)
        other.der.sort()
        other.dim += 1
        return other

    def string(self):
        """Returns a string representing the field
        """
        text = ""
        # derivatives acting on the field
        for mu in self.der:
            text += "d_" + str(mu) + " "
        # field itself
        text += self.name
        # spin indices
        if len(self.s1) + len(self.s2) > 0:
            text += "_"
            for a in self.s1:
                text += str(a)
            for b in self.s2:
                text += str(b)
        return text

    def latex(self):
        """Returns a latex expression representing the field
        """
        text = ""
        # derivatives acting on the field
        for mu in self.der:
            text += "\partial_" + str(mu) + " "
        # field itself
        text += self.name
        # spin indices
        if len(self.s1) + len(self.s2) > 0:
            text += "_{"
            for a in self.s1:
                text += str(a)
            for b in self.s2:
                text += str(b)
            text += "}"
        return text

    def __eq__(self, other):
        """Comparison operator
        """
        if (self.name == other.name
            and list(self.s1) == list(other.s1)
            and list(self.s2) == list(other.s2)
            and list(self.der) == list(other.der)
            and self.Grassmann == other.Grassmann
            and self.dim == other.dim):
            return True
        return False

    def __gt__(self, other):
        """Logical operator used for sorting the fields"""
        # (1) compare the names
        if self.name != other.name:
            return self.name < other.name
        # (2) if names are identical, compare the number of derivatives
        if len(self.der) != len(other.der):
            return len(self.der) > len(other.der)
        # (3) compare the indices of derivatives
        if self.der != other.der:
            return self.der > other.der
        # (4) compare the spin indices
        if self.s1 != other.s1:
            return self.s1 > other.s1
        return self.s2 > other.s2

##############################################################################

class SingleFieldChain:
    """Product of SingleField objects with a numerical factor
    """

    def __init__(self, name = None, s1 = [], s2 = [], dim = 0):
        """Constructor for a chain with zero or one terms
        If the name is fixed to 'None', then the chain is empty"""
        if name == None:
            self.field = [] # list of single fields in the chain
        else:
            self.field = [ SingleField(name, s1, s2, dim) ]
        self.factor = 1     # numerical (or symbolic) factor
        self.Grassmann = (len(s1) + len(s2)) % 2 == 1
                            # Grassmann number or not
        self.dim = dim      # scaling dimension

    def copy(self):
        """Creates a copy of the chain
        """
        other = SingleFieldChain(dim = self.dim)
        other.field = [field.copy() for field in self.field]
        other.factor = self.factor
        other.Grassmann = self.Grassmann
        return other

    def subs(self, x, y):
        """Returns a chain with x substituted by y
        """
        other = SingleFieldChain(dim = self.dim)
        other.field = [field.copy() for field in self.field]
        other.factor = self.factor.subs(x, y)
        other.Grassmann = self.Grassmann
        return other

    def check(self):
        """Consistency check of the scaling dimensions of the fields
                """
        dim, Grassmann = 0, False
        for field in self.field:
            dim += field.dim
            if field.Grassmann:
                Grassmann = not Grassmann
        if dim != self.dim:
            raise Error("SingleFieldChain::check()",
                        "Mismatch of the scaling dimension",
                        self.string())
        if Grassmann != self.Grassmann:
            raise Error("SingleFieldChain::check()",
                        "Mismatch of the spin statistics",
                        self.string())

    def string(self, sum = False):
        """Returns a string representing the chain of fields
        if the argument 'sum' is set to 'True', then
        the string is meant to be an element in a sum
        (not the first term) and is printed with a '+' in front"""
        text = ""
        c = self.factor
        if sum == True:
            if ((im(c) == 0 and re(c) < 0)
                or (re(c) == 0 and im(c) < 0)):
                text += " -"
                c = -c
            else:
                text += " +"
        if c != 1:
            text += " (" + str(c) + ")"
        for field in self.field:
            text += " " + field.string()
        return text

    def latex(self, sum = False):
        """Returns a string representing the chain of fields
        if the argument 'sum' is set to 'True', then
        the string is meant to be an element in a sum
        (not the first term) and is printed with a '+' in front"""
        text = ""
        c = self.factor
        if sum == True:
            if ((im(c) == 0 and re(c) < 0)
                or (re(c) == 0 and im(c) < 0)):
                text += " - "
                c = -c
            else:
                text += " +"
        if c != 1:
            text += " \\left(" + str(c) + "\\right)"
        for field in self.field:
            text += " " + field.latex()
        return text

    def sort(self):
        """Sorts the fields in the chain
        """
        minussign = False
        # use a simple bubble sort algorithm
        for i in range(len(self.field)):
            for j in range(len(self.field) - 1 - i):
                if self.field[j] > self.field[j + 1]:
                    self.field[j], self.field[j + 1] = (
                        self.field[j + 1], self.field[j]
                    )
                    if (self.field[j].Grassmann
                        and self.field[j + 1].Grassmann):
                        minussign = not minussign
                # the product of two identical Grassmann fields is zero
                elif (self.field[j] == self.field[j + 1]
                      and self.field[j].Grassmann):
                    self.factor = 0
                    return None
        if minussign:
            self.factor *= -1

    def d(self, mu):
        """Derivative acting on the term
        Note: this returns an array of SingleFieldChain
        """
        result = []
        n = len(self.field)
        for i in range(n):
            other = SingleFieldChain(dim = self.dim + 1)
            other.field = (
                [self.field[j].copy() for j in range(i)]
                + [ self.field[i].d(mu) ]
                + [self.field[j].copy() for j in range(i+1, n)]
            )
            other.factor = self.factor
            other.Grassmann = self.Grassmann
            result.append(other)
        return result

    def __rmul__(self, c):
        """Multiplication by a real number
        """
        result = self.copy()
        result.factor *= c
        return result

    def __mul__(self, other):
        """Multiplication of two chains of fields
        """
        result = SingleFieldChain(dim = self.dim + other.dim)
        result.field = self.field + other.field
        result.factor = self.factor * other.factor
        result.Grassmann = self.Grassmann != other.Grassmann
        return result

##############################################################################

class Field:
    """Sum of SingleFieldChain objects
    """

    def __init__(self, name = None, s1 = [], s2 = [], dim = 0,
                 Grassmann = None):
        """Constructor for a single field or an empty expression
        """
        if name == None:
            self.term = []  # list of chains in the sum
            if Grassmann == None:
                raise Error("Field::__init__",
                            "The spin statistics of an empty "
                            "field must be specified", "")
            self.Grassmann = Grassmann
        else:
            self.term = [ SingleFieldChain(name, s1, s2, dim) ]
            if Grassmann != None:
                raise Error("Field::__init__",
                            "The spin statistics of a non-empty "
                            "field cannot be specified",
                            self.string())
            self.Grassmann = (len(s1) + len(s2)) % 2 == 1
                            # Grassmann number or not
        self.dim = dim      # scaling dimension

    def copy(self):
        """Creates a copy of the field
        """
        other = Field(dim = self.dim,
                      Grassmann = self.Grassmann)
        other.term = [term.copy() for term in self.term]
        return other

    def subs(self, x, y):
        """Returns a field with x substituted by y
        """
        other = Field(dim = self.dim,
                      Grassmann = self.Grassmann)
        other.term = [term.subs(x, y) for term in self.term]
        return other

    def check(self):
        """Consistency check of the scaling dimensions of the fields
        """
        for term in self.term:
            term.check()
            if term.dim != self.dim:
                raise Error("Field::check()",
                            "Mismatch of the scaling dimension",
                            self.string())
            if term.Grassmann != self.Grassmann:
                raise Error("Field::check()",
                            "Mismatch of the spin statistics",
                            self.string())

    def iszero(self):
        """Check if a field is zero
        """
        return len(self.term) == 0

    def string(self):
        """Returns a string representing the field
        """
        if len(self.term) == 0:
            return " 0"
        text = self.term[0].string()
        for i in range(1, len(self.term)):
            text += self.term[i].string(sum = True)
        return text

    def latex(self):
        """Returns a latex expression representing the field
        """
        if len(self.term) == 0:
            return "0"
        text = self.term[0].latex()
        for i in range(1, len(self.term)):
            text += self.term[i].latex(sum = True)
        return text

    def print(self, latexoutput = latexoutput):
        """Print the field in the chosen format"
        """
        if latexoutput:
            print(self.latex())
        else:
            print(self.string())

    def sort(self):
        """Sorts the various terms
        """
        # (1) sort the fields in each term
        for term in self.term:
            term.sort()
        # (2) sort the terms in the sum
        self.term.sort(key = lambda term: term.field)
        # (3) combine identical terms
        for i in reversed(range(1,len(self.term))):
            if self.term[i-1].field == self.term[i].field:
                self.term[i-1].factor += self.term[i].factor
                self.term.pop(i)
        # (4) simplify the factors and remove zeroes
        for i in reversed(range(len(self.term))):
            self.term[i].factor = simplify(self.term[i].factor)
            if self.term[i].factor == 0:
                self.term.pop(i)

    def d(self, mu):
        """Derivative acting on the expression
        """
        other = Field(dim = self.dim + 1,
                      Grassmann = self.Grassmann)
        for term in self.term:
            other.term += term.d(mu)
        return other
        
    def dd(self, a, b):
        """Alternate definition of the derivative using spin indices
        """
        return sum(
            lambda aa, bb, mu:
                sigma[mu[0]][a][b],
            lambda aa, bb, mu:
                self.d(mu[0]),
            0, 0, 1
        )

    def box(self):
        """Laplacian of the field
        """
        return (self.d(0).d(0) - self.d(1).d(1)
                - self.d(2).d(2) - self.d(3).d(3))

    def __add__(self, other):
        """Sum of two expressions
        """
        if self.dim != other.dim:
            raise Error("Field::__add__(other)",
                        "Adding terms with different scaling dimensions",
                        self.string())
        if self.Grassmann != other.Grassmann:
            raise Error("Field::__add__(other)",
                        "Adding terms with different spin statistics",
                        self.string())
        result = Field(dim = self.dim,
                       Grassmann = self.Grassmann)
        result.term = ([term.copy() for term in self.term]
            + [term.copy() for term in other.term])
        return result

    def __rmul__(self, c):
        """Used to define multiplication by a real number
        """
        other = Field(dim = self.dim,
                      Grassmann = self.Grassmann)
        other.term = [c * term.copy() for term in self.term]
        return other

    def __sub__(self, other):
        """Difference of two expressions
        """
        return self + (-1) * other

    def __mul__(self, other):
        """Multiplication of two fields
        """
        result = Field(dim = self.dim + other.dim,
                       Grassmann = self.Grassmann != other.Grassmann)
        result.term = [term1 * term2
                       for term1 in self.term
                       for term2 in other.term]
        return result

    def coefficientset(self):
        """Set of coefficients of the terms"""
        coeffs = set()
        for term in self.term:
            coeffs.add(term.factor)
        return coeffs

    def bosonicpart(self):
        """Returns the field stripped from all its fermions
        """
        result = Field(dim = self.dim, Grassmann = self.Grassmann)
        for term in self.term:
            flag = True
            for field in term.field:
                if field.Grassmann == 1:
                    flag = False
#                if (field.name == "f"
#                    or field.name == "f^*"
#                    or field.name == "theta"):
#                    flag = False
            if flag:
                result.term.append(term.copy())
        return result

    def staticpart(self):
        """Returns the field stripped from all its time derivatives
        """
        result = Field(dim = self.dim, Grassmann = self.Grassmann)
        for term in self.term:
            flag = True
            for field in term.field:
                if 0 in field.der:
                    flag = False
                if field.name == "A_0":
                    flag = False
            if flag:
                result.term.append(term.copy())
        return result

    def WessZuminogauge(self):
        """Returns the field in the Wess-Zumino gauge,
        that is with c = chi = N = 0
        """
        result = Field(dim = self.dim, Grassmann = self.Grassmann)
        for term in self.term:
            flag = True
            for field in term.field:
                if (field.name == "c"
                    or field.name == "\\chi"
                    or field.name == "\\chi^*"
                    or field.name == "N"
                    or field.name == "N^*"):
                    flag = False
            if flag:
                result.term.append(term.copy())
        return result

##############################################################################

# tables used for checking the scaling dimension and
# spin statistics of the supermultiplet
dimTable = ( (0, half, half, 1),
             (half, 1, 1, 1+half),
             (half, 1, 1, 1+half),
             (1, 1+half, 1+half, 2) )
GrassmannTable = ( (False, True, True, False),
                   (True, False, False, True),
                   (True, False, False, True),
                   (False, True, True, False) )

# rules for superfield multiplication and covariant derivatives
MultRule = ((-1, 0, 3, 0, 1, 0, 2),
             (-1, 1, 1, 1, 0, 0, 1),
             (-1, 1, 2, 1, 0, 0, 2),
             (-1, 1, 3, 0, 1, 1, 2),
             (-1, 1, 3, 1, 2, 0, 1),
             (-1, 2, 1, 2, 0, 0, 1),
             (-1, 2, 2, 2, 0, 0, 2),
             (-1, 2, 3, 0, 1, 2, 2),
             (-1, 2, 3, 2, 2, 0, 1),
             (-1, 3, 0, 1, 0, 2, 0),
             (-1, 3, 1, 1, 1, 2, 0),
             (-1, 3, 1, 2, 0, 1, 1),
             (-1, 3, 2, 1, 2, 2, 0),
             (-1, 3, 2, 2, 0, 1, 2),
             (-1, 3, 3, 0, 1, 3, 2),
             (-1, 3, 3, 1, 0, 2, 3),
             (-1, 3, 3, 1, 1, 2, 2),
             (-1, 3, 3, 1, 3, 2, 0),
             (-1, 3, 3, 2, 2, 1, 1),
             (-1, 3, 3, 3, 1, 0, 2),
             (1, 0, 0, 0, 0, 0, 0),
             (1, 0, 1, 0, 0, 0, 1),
             (1, 0, 1, 0, 1, 0, 0),
             (1, 0, 2, 0, 0, 0, 2),
             (1, 0, 2, 0, 2, 0, 0),
             (1, 0, 3, 0, 0, 0, 3),
             (1, 0, 3, 0, 2, 0, 1),
             (1, 0, 3, 0, 3, 0, 0),
             (1, 1, 0, 0, 0, 1, 0),
             (1, 1, 0, 1, 0, 0, 0),
             (1, 1, 1, 0, 0, 1, 1),
             (1, 1, 1, 0, 1, 1, 0),
             (1, 1, 1, 1, 1, 0, 0),
             (1, 1, 2, 0, 0, 1, 2),
             (1, 1, 2, 0, 2, 1, 0),
             (1, 1, 2, 1, 2, 0, 0),
             (1, 1, 3, 0, 0, 1, 3),
             (1, 1, 3, 0, 2, 1, 1),
             (1, 1, 3, 0, 3, 1, 0),
             (1, 1, 3, 1, 0, 0, 3),
             (1, 1, 3, 1, 1, 0, 2),
             (1, 1, 3, 1, 3, 0, 0),
             (1, 2, 0, 0, 0, 2, 0),
             (1, 2, 0, 2, 0, 0, 0),
             (1, 2, 1, 0, 0, 2, 1),
             (1, 2, 1, 0, 1, 2, 0),
             (1, 2, 1, 2, 1, 0, 0),
             (1, 2, 2, 0, 0, 2, 2),
             (1, 2, 2, 0, 2, 2, 0),
             (1, 2, 2, 2, 2, 0, 0),
             (1, 2, 3, 0, 0, 2, 3),
             (1, 2, 3, 0, 2, 2, 1),
             (1, 2, 3, 0, 3, 2, 0),
             (1, 2, 3, 2, 0, 0, 3),
             (1, 2, 3, 2, 1, 0, 2),
             (1, 2, 3, 2, 3, 0, 0),
             (1, 3, 0, 0, 0, 3, 0),
             (1, 3, 0, 2, 0, 1, 0),
             (1, 3, 0, 3, 0, 0, 0),
             (1, 3, 1, 0, 0, 3, 1),
             (1, 3, 1, 0, 1, 3, 0),
             (1, 3, 1, 1, 0, 2, 1),
             (1, 3, 1, 2, 1, 1, 0),
             (1, 3, 1, 3, 0, 0, 1),
             (1, 3, 1, 3, 1, 0, 0),
             (1, 3, 2, 0, 0, 3, 2),
             (1, 3, 2, 0, 2, 3, 0),
             (1, 3, 2, 1, 0, 2, 2),
             (1, 3, 2, 2, 2, 1, 0),
             (1, 3, 2, 3, 0, 0, 2),
             (1, 3, 2, 3, 2, 0, 0),
             (1, 3, 3, 0, 0, 3, 3),
             (1, 3, 3, 0, 2, 3, 1),
             (1, 3, 3, 0, 3, 3, 0),
             (1, 3, 3, 1, 2, 2, 1),
             (1, 3, 3, 2, 0, 1, 3),
             (1, 3, 3, 2, 1, 1, 2),
             (1, 3, 3, 2, 3, 1, 0),
             (1, 3, 3, 3, 0, 0, 3),
             (1, 3, 3, 3, 2, 0, 1),
             (1, 3, 3, 3, 3, 0, 0))
MultGrassmannRule = ((-1, 0, 1, 0, 1, 0, 0),
                     (-1, 0, 2, 0, 2, 0, 0),
                     (-1, 0, 3, 0, 2, 0, 1),
                     (-1, 1, 0, 1, 0, 0, 0),
                     (-1, 1, 1, 0, 1, 1, 0),
                     (-1, 1, 2, 0, 2, 1, 0),
                     (-1, 1, 3, 0, 2, 1, 1),
                     (-1, 1, 3, 1, 0, 0, 3),
                     (-1, 1, 3, 1, 2, 0, 1),
                     (-1, 1, 3, 1, 3, 0, 0),
                     (-1, 2, 0, 2, 0, 0, 0),
                     (-1, 2, 1, 0, 1, 2, 0),
                     (-1, 2, 2, 0, 2, 2, 0),
                     (-1, 2, 3, 0, 2, 2, 1),
                     (-1, 2, 3, 2, 0, 0, 3),
                     (-1, 2, 3, 2, 2, 0, 1),
                     (-1, 2, 3, 2, 3, 0, 0),
                     (-1, 3, 0, 2, 0, 1, 0),
                     (-1, 3, 1, 0, 1, 3, 0),
                     (-1, 3, 1, 1, 0, 2, 1),
                     (-1, 3, 1, 1, 1, 2, 0),
                     (-1, 3, 1, 3, 1, 0, 0),
                     (-1, 3, 2, 0, 2, 3, 0),
                     (-1, 3, 2, 1, 0, 2, 2),
                     (-1, 3, 2, 1, 2, 2, 0),
                     (-1, 3, 2, 3, 2, 0, 0),
                     (-1, 3, 3, 0, 2, 3, 1),
                     (-1, 3, 3, 1, 1, 2, 2),
                     (-1, 3, 3, 2, 0, 1, 3),
                     (-1, 3, 3, 2, 2, 1, 1),
                     (-1, 3, 3, 2, 3, 1, 0),
                     (-1, 3, 3, 3, 2, 0, 1),
                     (1, 0, 0, 0, 0, 0, 0),
                     (1, 0, 1, 0, 0, 0, 1),
                     (1, 0, 2, 0, 0, 0, 2),
                     (1, 0, 3, 0, 0, 0, 3),
                     (1, 0, 3, 0, 1, 0, 2),
                     (1, 0, 3, 0, 3, 0, 0),
                     (1, 1, 0, 0, 0, 1, 0),
                     (1, 1, 1, 0, 0, 1, 1),
                     (1, 1, 1, 1, 0, 0, 1),
                     (1, 1, 1, 1, 1, 0, 0),
                     (1, 1, 2, 0, 0, 1, 2),
                     (1, 1, 2, 1, 0, 0, 2),
                     (1, 1, 2, 1, 2, 0, 0),
                     (1, 1, 3, 0, 0, 1, 3),
                     (1, 1, 3, 0, 1, 1, 2),
                     (1, 1, 3, 0, 3, 1, 0),
                     (1, 1, 3, 1, 1, 0, 2),
                     (1, 2, 0, 0, 0, 2, 0),
                     (1, 2, 1, 0, 0, 2, 1),
                     (1, 2, 1, 2, 0, 0, 1),
                     (1, 2, 1, 2, 1, 0, 0),
                     (1, 2, 2, 0, 0, 2, 2),
                     (1, 2, 2, 2, 0, 0, 2),
                     (1, 2, 2, 2, 2, 0, 0),
                     (1, 2, 3, 0, 0, 2, 3),
                     (1, 2, 3, 0, 1, 2, 2),
                     (1, 2, 3, 0, 3, 2, 0),
                     (1, 2, 3, 2, 1, 0, 2),
                     (1, 3, 0, 0, 0, 3, 0),
                     (1, 3, 0, 1, 0, 2, 0),
                     (1, 3, 0, 3, 0, 0, 0),
                     (1, 3, 1, 0, 0, 3, 1),
                     (1, 3, 1, 2, 0, 1, 1),
                     (1, 3, 1, 2, 1, 1, 0),
                     (1, 3, 1, 3, 0, 0, 1),
                     (1, 3, 2, 0, 0, 3, 2),
                     (1, 3, 2, 2, 0, 1, 2),
                     (1, 3, 2, 2, 2, 1, 0),
                     (1, 3, 2, 3, 0, 0, 2),
                     (1, 3, 3, 0, 0, 3, 3),
                     (1, 3, 3, 0, 1, 3, 2),
                     (1, 3, 3, 0, 3, 3, 0),
                     (1, 3, 3, 1, 0, 2, 3),
                     (1, 3, 3, 1, 2, 2, 1),
                     (1, 3, 3, 1, 3, 2, 0),
                     (1, 3, 3, 2, 1, 1, 2),
                     (1, 3, 3, 3, 0, 0, 3),
                     (1, 3, 3, 3, 1, 0, 2),
                     (1, 3, 3, 3, 3, 0, 0))
covDRule = (((-1, 0, 0, 0, 1, 0),
             (-1, 0, 3, 0, 1, 3),
             (-1, 2, 1, 0, 3, 1),
             (-1, 2, 2, 0, 3, 2),
             (-I, 0, 3, 1, 0, 2),
             (-I, 0, 3, 4, 0, 2),
             (-I, 1, 3, 1, 1, 2),
             (-I, 1, 3, 4, 1, 2),
             (-I, 2, 3, 1, 2, 2),
             (-I, 2, 3, 4, 2, 2),
             (-I, 3, 3, 1, 3, 2),
             (-I, 3, 3, 4, 3, 2),
             (I, 0, 1, 1, 0, 0),
             (I, 0, 1, 4, 0, 0),
             (I, 0, 2, 2, 0, 0),
             (I, 0, 3, 2, 0, 1),
             (I, 1, 1, 1, 1, 0),
             (I, 1, 1, 4, 1, 0),
             (I, 1, 2, 2, 1, 0),
             (I, 1, 3, 2, 1, 1),
             (I, 2, 1, 1, 2, 0),
             (I, 2, 1, 4, 2, 0),
             (I, 2, 2, 2, 2, 0),
             (I, 2, 3, 2, 2, 1),
             (I, 3, 1, 1, 3, 0),
             (I, 3, 1, 4, 3, 0),
             (I, 3, 2, 2, 3, 0),
             (I, 3, 3, 2, 3, 1),
             (1, 0, 1, 0, 1, 1),
             (1, 0, 2, 0, 1, 2),
             (1, 0, 2, 3, 0, 0),
             (1, 0, 3, 3, 0, 1),
             (1, 1, 2, 3, 1, 0),
             (1, 1, 3, 3, 1, 1),
             (1, 2, 0, 0, 3, 0),
             (1, 2, 2, 3, 2, 0),
             (1, 2, 3, 0, 3, 3),
             (1, 2, 3, 3, 2, 1),
             (1, 3, 2, 3, 3, 0),
             (1, 3, 3, 3, 3, 1)),
            ((-1, 0, 0, 0, 2, 0),
             (-1, 0, 1, 3, 0, 0),
             (-1, 0, 3, 0, 2, 3),
             (-1, 1, 0, 0, 3, 0),
             (-1, 1, 1, 3, 1, 0),
             (-1, 1, 3, 0, 3, 3),
             (-1, 2, 1, 3, 2, 0),
             (-1, 3, 1, 3, 3, 0),
             (-I, 0, 2, 4, 0, 0),
             (-I, 0, 3, 2, 0, 2),
             (-I, 0, 3, 4, 0, 1),
             (-I, 1, 2, 4, 1, 0),
             (-I, 1, 3, 2, 1, 2),
             (-I, 1, 3, 4, 1, 1),
             (-I, 2, 2, 4, 2, 0),
             (-I, 2, 3, 2, 2, 2),
             (-I, 2, 3, 4, 2, 1),
             (-I, 3, 2, 4, 3, 0),
             (-I, 3, 3, 2, 3, 2),
             (-I, 3, 3, 4, 3, 1),
             (I, 0, 1, 2, 0, 0),
             (I, 0, 2, 1, 0, 0),
             (I, 0, 3, 1, 0, 1),
             (I, 1, 1, 2, 1, 0),
             (I, 1, 2, 1, 1, 0),
             (I, 1, 3, 1, 1, 1),
             (I, 2, 1, 2, 2, 0),
             (I, 2, 2, 1, 2, 0),
             (I, 2, 3, 1, 2, 1),
             (I, 3, 1, 2, 3, 0),
             (I, 3, 2, 1, 3, 0),
             (I, 3, 3, 1, 3, 1),
             (1, 0, 1, 0, 2, 1),
             (1, 0, 2, 0, 2, 2),
             (1, 0, 3, 3, 0, 2),
             (1, 1, 1, 0, 3, 1),
             (1, 1, 2, 0, 3, 2),
             (1, 1, 3, 3, 1, 2),
             (1, 2, 3, 3, 2, 2),
             (1, 3, 3, 3, 3, 2)))
covDbarRule = (((-1, 0, 2, 0, 0, 3),
                 (-1, 1, 2, 0, 1, 3),
                 (-1, 2, 1, 3, 0, 1),
                 (-1, 2, 2, 0, 2, 3),
                 (-1, 2, 2, 3, 0, 2),
                 (-1, 3, 1, 3, 1, 1),
                 (-1, 3, 2, 0, 3, 3),
                 (-1, 3, 2, 3, 1, 2),
                 (-I, 1, 0, 1, 0, 0),
                 (-I, 1, 0, 4, 0, 0),
                 (-I, 1, 3, 1, 0, 3),
                 (-I, 1, 3, 4, 0, 3),
                 (-I, 2, 0, 2, 0, 0),
                 (-I, 2, 3, 2, 0, 3),
                 (-I, 3, 0, 2, 1, 0),
                 (-I, 3, 1, 1, 2, 1),
                 (-I, 3, 1, 4, 2, 1),
                 (-I, 3, 2, 1, 2, 2),
                 (-I, 3, 2, 4, 2, 2),
                 (-I, 3, 3, 2, 1, 3),
                 (I, 1, 1, 1, 0, 1),
                 (I, 1, 1, 4, 0, 1),
                 (I, 1, 2, 1, 0, 2),
                 (I, 1, 2, 4, 0, 2),
                 (I, 2, 1, 2, 0, 1),
                 (I, 2, 2, 2, 0, 2),
                 (I, 3, 0, 1, 2, 0),
                 (I, 3, 0, 4, 2, 0),
                 (I, 3, 1, 2, 1, 1),
                 (I, 3, 2, 2, 1, 2),
                 (I, 3, 3, 1, 2, 3),
                 (I, 3, 3, 4, 2, 3),
                 (1, 0, 0, 0, 0, 1),
                 (1, 1, 0, 0, 1, 1),
                 (1, 2, 0, 0, 2, 1),
                 (1, 2, 0, 3, 0, 0),
                 (1, 2, 3, 3, 0, 3),
                 (1, 3, 0, 0, 3, 1),
                 (1, 3, 0, 3, 1, 0),
                 (1, 3, 3, 3, 1, 3)),
                ((-1, 1, 0, 3, 0, 0),
                 (-1, 1, 3, 3, 0, 3),
                 (-1, 3, 1, 3, 2, 1),
                 (-1, 3, 2, 3, 2, 2),
                 (-I, 1, 0, 2, 0, 0),
                 (-I, 1, 3, 2, 0, 3),
                 (-I, 2, 0, 1, 0, 0),
                 (-I, 2, 1, 4, 0, 1),
                 (-I, 2, 2, 4, 0, 2),
                 (-I, 2, 3, 1, 0, 3),
                 (-I, 3, 0, 1, 1, 0),
                 (-I, 3, 1, 2, 2, 1),
                 (-I, 3, 1, 4, 1, 1),
                 (-I, 3, 2, 2, 2, 2),
                 (-I, 3, 2, 4, 1, 2),
                 (-I, 3, 3, 1, 1, 3),
                 (I, 1, 1, 2, 0, 1),
                 (I, 1, 2, 2, 0, 2),
                 (I, 2, 0, 4, 0, 0),
                 (I, 2, 1, 1, 0, 1),
                 (I, 2, 2, 1, 0, 2),
                 (I, 2, 3, 4, 0, 3),
                 (I, 3, 0, 2, 2, 0),
                 (I, 3, 0, 4, 1, 0),
                 (I, 3, 1, 1, 1, 1),
                 (I, 3, 2, 1, 1, 2),
                 (I, 3, 3, 2, 2, 3),
                 (I, 3, 3, 4, 1, 3),
                 (1, 0, 0, 0, 0, 2),
                 (1, 0, 1, 0, 0, 3),
                 (1, 1, 0, 0, 1, 2),
                 (1, 1, 1, 0, 1, 3),
                 (1, 1, 1, 3, 0, 1),
                 (1, 1, 2, 3, 0, 2),
                 (1, 2, 0, 0, 2, 2),
                 (1, 2, 1, 0, 2, 3),
                 (1, 3, 0, 0, 3, 2),
                 (1, 3, 0, 3, 2, 0),
                 (1, 3, 1, 0, 3, 3),
                 (1, 3, 3, 3, 2, 3)))

class SuperField:
    """SUSY multiplet"""

    def __init__(self, type=None, dim=1, Grassmann=None):
        """Constructor: superfield of type
        chiral, anti-chiral, vector, or generic superfield with index
        """
        if (type == "chiral"): # chiral supermultiplet
            self.Grassmann = False
            self.dim = dim  # scaling dimension
            phi = Field("\\phi", dim = dim)
            chi = [Field("\\chi", [a], [], dim = dim + half)
                   for a in range(1, 3)]
            F = Field("F", dim = dim + 1)
            self.component = [[phi,
                               Field(dim = dim + half,
                                     Grassmann = True),
                               Field(dim = dim + half,
                                     Grassmann = True),
                               Field(dim = dim + 1,
                                     Grassmann = False)],
                              [-2 * chi[0],
                               I * phi.d(0) + I * phi.d(3),
                               I * phi.d(1) + phi.d(2),
                               Field(dim = dim + 1 + half,
                                     Grassmann = True)],
                              [-2 * chi[1],
                               I * phi.d(1) - phi.d(2),
                               I * phi.d(0) - I * phi.d(3),
                               Field(dim = dim + 1 + half,
                                     Grassmann = True)],
                              [-2 * F,
                               2 * (I * chi[1].d(0) - I * chi[0].d(1)
                               + chi[0].d(2) + I * chi[1].d(3)),
                               2 * (-I * chi[0].d(0) + I * chi[1].d(1)
                               + chi[1].d(2) + I * chi[0].d(3)),
                               phi.box() ]]
        elif (type == "antichiral"): # antichiral supermultiplet
            self.Grassmann = False
            self.dim = dim  # scaling dimension
            phiC = Field("\\phi^*", dim = dim)
            chiC = [Field("\\chi^*", [], [b], dim = dim + half)
                   for b in range(1, 3)]
            FC = Field("F^*", dim = dim + 1)
            self.component = [[phiC,
                               2 * chiC[0],
                               2 * chiC[1],
                               2 * FC],
                              [Field(dim = dim + half,
                                     Grassmann = True),
                               -I * phiC.d(0) - I * phiC.d(3),
                               -I * phiC.d(1) - phiC.d(2),
                               2 * (-I * chiC[1].d(0) + I * chiC[0].d(1)
                                + chiC[0].d(2) - I * chiC[1].d(3))],
                              [Field(dim = dim + half,
                                     Grassmann = True),
                               -I * phiC.d(1) + phiC.d(2),
                               -I * phiC.d(0) + I * phiC.d(3),
                               2 * (I * chiC[0].d(0) - I * chiC[1].d(1)
                                + chiC[1].d(2) - I * chiC[0].d(3))
                               ],
                              [Field(dim = dim + 1,
                                     Grassmann = False),
                               Field(dim=dim + 1 + half,
                                     Grassmann=True),
                               Field(dim=dim + 1 + half,
                                     Grassmann=True),
                               phiC.box() ]]
        elif (type == "vector"): # vector supermultiplet
            self.Grassmann = False
            self.dim = dim - 1  # scaling dimension
            c = Field("c", dim = dim - 1)
            chi = [Field("\\chi", [a], [], dim=dim - half)
                    for a in range(1, 3)]
            chiC = [Field("\\chi^*", [], [b], dim=dim - half)
                    for b in range(1, 3)]
            N = Field("N", dim = dim)
            NC = Field("N^*", dim = dim)
            Amu = [Field("A_" + str(mu), [], [], dim = dim)
                   for mu in range(4)]
            A = [[sum(lambda aa, bb, mu:
                        half * eta[mu[0]][mu[1]] * sigma[mu[1]][a][b],
                      lambda aa, bb, mu: Amu[mu[0]],
                      0, 0, 2)
                  for b in range(2)]
                 for a in range(2)]
            # A = [[Field("A", [a], [b], dim = dim)
            #       for b in range(1,3)]
            #      for a in range(1,3)]
            lam = [Field("\\lambda", [a], [], dim=dim + half)
                   for a in range(1, 3)]
            lamC = [Field("\\lambda^*", [], [b], dim=dim + half)
                    for b in range(1, 3)]
            D = Field("D", dim = dim + 1)
            self.component = [[c,
                               chiC[0],
                               chiC[1],
                               2 * NC ],
                              [-1 * chi[0],
                               A[0][0],
                               A[0][1],
                               (lam[0]
                                - I * chiC[1].d(0)
                                + I * chiC[0].d(1)
                                + chiC[0].d(2)
                                - I * chiC[1].d(3)) ],
                              [-1 * chi[1],
                               A[1][0],
                               A[1][1],
                               (lam[1]
                                + I * chiC[0].d(0)
                                - I * chiC[1].d(1)
                                + chiC[1].d(2)
                                - I * chiC[0].d(3)) ],
                              [-2 * N,
                               (lamC[0]
                                + I * chi[1].d(0)
                                - I * chi[0].d(1)
                                + chi[0].d(2)
                                + I * chi[1].d(3)),
                               (lamC[1]
                                - I * chi[0].d(0)
                                + I * chi[1].d(1)
                                + chi[1].d(2)
                                + I * chi[0].d(3)),
                               (-1) * D + c.box()]]
        elif (type == "WZvector"): # vector supermultiplet in the Wess-Zumino gauge
            self.Grassmann = False
            self.dim = dim - 1  # scaling dimension
            Amu = [Field("A_" + str(mu), [], [], dim = dim)
                   for mu in range(4)]
            A = [[sum(lambda aa, bb, mu:
                        half * eta[mu[0]][mu[1]] * sigma[mu[1]][a][b],
                      lambda aa, bb, mu: Amu[mu[0]],
                      0, 0, 2)
                  for b in range(2)]
                 for a in range(2)]
            lam = [Field("\\lambda", [a], [], dim=dim + half)
                   for a in range(1, 3)]
            lamC = [Field("\\lambda^*", [], [b], dim=dim + half)
                    for b in range(1, 3)]
            D = Field("D", dim = dim + 1)
            self.component = [[Field(dim=dim - 1, Grassmann=False),
                               Field(dim=dim - half, Grassmann=True),
                               Field(dim=dim - half, Grassmann=True),
                               Field(dim=dim, Grassmann=False) ],
                              [Field(dim=dim - half, Grassmann=True),
                               A[0][0],
                               A[0][1],
                               lam[0]],
                              [Field(dim=dim - half, Grassmann=True),
                               A[1][0],
                               A[1][1],
                               lam[1]],
                              [Field(dim=dim, Grassmann=False),
                               lamC[0],
                               lamC[1],
                               (-1) * D]]
        elif (type == "chiralcoupling"): # chiral supermultiplet
            self.Grassmann = False
            self.dim = 0  # scaling dimension
            g = Field("g", dim = 0)
            theta = Field("theta", dim = 0)
            phi = half * sqrt(2) * (g + I * theta)
            chi = [Field("\\xi", [a], [], dim = half)
                   for a in range(1, 3)]
            F = Field("f", dim = 1)
            self.component = [[phi,
                               Field(dim = half,
                                     Grassmann = True),
                               Field(dim = half,
                                     Grassmann = True),
                               Field(dim = 1,
                                     Grassmann = False)],
                              [-2 * chi[0],
                               I * phi.d(0) + I * phi.d(3),
                               I * phi.d(1) + phi.d(2),
                               Field(dim = 1 + half,
                                     Grassmann = True)],
                              [-2 * chi[1],
                               I * phi.d(1) - phi.d(2),
                               I * phi.d(0) - I * phi.d(3),
                               Field(dim = 1 + half,
                                     Grassmann = True)],
                              [-2 * F,
                               2 * (I * chi[1].d(0) - I * chi[0].d(1)
                               + chi[0].d(2) + I * chi[1].d(3)),
                               2 * (-I * chi[0].d(0) + I * chi[1].d(1)
                               + chi[1].d(2) + I * chi[0].d(3)),
                               phi.box() ]]
        elif (type == "antichiralcoupling"): # antichiral supermultiplet
            self.Grassmann = False
            self.dim = 0  # scaling dimension
            g = Field("g", dim = 0)
            theta = Field("theta", dim = 0)
            phiC = half * sqrt(2) * (g - I * theta)
            chiC = [Field("\\xi^*", [], [b], dim = half)
                   for b in range(1, 3)]
            FC = Field("f^*", dim = 1)
            self.component = [[phiC,
                               2 * chiC[0],
                               2 * chiC[1],
                               2 * FC],
                              [Field(dim = half,
                                     Grassmann = True),
                               -I * phiC.d(0) - I * phiC.d(3),
                               -I * phiC.d(1) - phiC.d(2),
                               2 * (-I * chiC[1].d(0) + I * chiC[0].d(1)
                                + chiC[0].d(2) - I * chiC[1].d(3))],
                              [Field(dim = half,
                                     Grassmann = True),
                               -I * phiC.d(1) + phiC.d(2),
                               -I * phiC.d(0) + I * phiC.d(3),
                               2 * (I * chiC[0].d(0) - I * chiC[1].d(1)
                                + chiC[1].d(2) - I * chiC[0].d(3))
                               ],
                              [Field(dim = 1,
                                     Grassmann = False),
                               Field(dim=1 + half,
                                     Grassmann=True),
                               Field(dim=1 + half,
                                     Grassmann=True),
                               phiC.box() ]]
        # elif (type == "realscalar"): # real scalar field put in a multiplet
        #     self.Grassmann = False
        #     self.dim = dim  # scaling dimension
        #     g = Field("g", dim = dim)
        #     self.component = [[g,
        #                        Field(dim = dim + half,
        #                              Grassmann = True),
        #                        Field(dim = dim + half,
        #                              Grassmann = True),
        #                        Field(dim = dim + 1,
        #                              Grassmann = False)],
        #                       [Field(dim = dim + half,
        #                              Grassmann = True),
        #                        Field(dim=dim + 1,
        #                              Grassmann=False),
        #                        Field(dim=dim + 1,
        #                              Grassmann=False),
        #                        Field(dim = dim + 1 + half,
        #                              Grassmann = True)],
        #                       [Field(dim = dim + half,
        #                              Grassmann = True),
        #                        Field(dim=dim + 1,
        #                              Grassmann=False),
        #                        Field(dim=dim + 1,
        #                              Grassmann=False),
        #                        Field(dim = dim + 1 + half,
        #                              Grassmann = True)],
        #                       [Field(dim = dim + 1,
        #                              Grassmann = False),
        #                        Field(dim=dim + 1 + half,
        #                              Grassmann=True),
        #                        Field(dim=dim + 1 + half,
        #                              Grassmann=True),
        #                        g.box() ]]
        # elif (type == "chiralscalar"): # chiral supermultiplet with only the scalar part non-zero
        #     self.Grassmann = False
        #     self.dim = dim  # scaling dimension
        #     g = Field("g", dim = dim)
        #     self.component = [[g,
        #                        Field(dim = dim + half,
        #                              Grassmann = True),
        #                        Field(dim = dim + half,
        #                              Grassmann = True),
        #                        Field(dim = dim + 1,
        #                              Grassmann = False)],
        #                       [Field(dim = dim + half,
        #                              Grassmann = True),
        #                        I * g.d(0) + I * g.d(3),
        #                        I * g.d(1) + g.d(2),
        #                        Field(dim = dim + 1 + half,
        #                              Grassmann = True)],
        #                       [Field(dim = dim + half,
        #                              Grassmann = True),
        #                        I * g.d(1) - g.d(2),
        #                        I * g.d(0) - I * g.d(3),
        #                        Field(dim = dim + 1 + half,
        #                              Grassmann = True)],
        #                       [Field(dim = dim + 1,
        #                              Grassmann = False),
        #                        Field(dim=dim + 1 + half,
        #                              Grassmann=True),
        #                        Field(dim=dim + 1 + half,
        #                              Grassmann=True),
        #                        g.box() ]]
        # elif (type == "antichiralscalar"): # antichiral supermultiplet with only the scalar part non-zero
        #     self.Grassmann = False
        #     self.dim = dim  # scaling dimension
        #     g = Field("g", dim = dim)
        #     self.component = [[g,
        #                        Field(dim=dim + half,
        #                              Grassmann=True),
        #                        Field(dim=dim + half,
        #                              Grassmann=True),
        #                        Field(dim=dim + 1,
        #                              Grassmann=False)],
        #                       [Field(dim = dim + half,
        #                              Grassmann = True),
        #                        -I * g.d(0) - I * g.d(3),
        #                        -I * g.d(1) - g.d(2),
        #                        Field(dim=dim + 1 + half,
        #                              Grassmann=True)],
        #                       [Field(dim = dim + half,
        #                              Grassmann = True),
        #                        -I * g.d(1) + g.d(2),
        #                        -I * g.d(0) + I * g.d(3),
        #                        Field(dim=dim + 1 + half,
        #                              Grassmann=True)],
        #                       [Field(dim = dim + 1,
        #                              Grassmann = False),
        #                        Field(dim=dim + 1 + half,
        #                              Grassmann=True),
        #                        Field(dim=dim + 1 + half,
        #                              Grassmann=True),
        #                        g.box() ]]
        elif (type == None):
            if Grassmann == None:
                raise Error("SuperField::__init__",
                            "The spin statistics of an empty "
                            "superfield must be specified", "")
            self.Grassmann = Grassmann
            self.dim = dim  # scaling dimension
            self.component = [ [ Field(dim = dimTable[i][j] + dim,
                                       Grassmann = (GrassmannTable[i][j]
                                                    != Grassmann))
                                 for j in range(4) ]
                               for i in range(4) ]
        else:  # generic supermultiplet with index "type"
            self.Grassmann = False
            self.dim = dim  # scaling dimension
            self.component = [[Field("a" + str(type),
                                     dim = dim),
                               Field("\\chi" + str(type),
                                     [], [1], dim + half),
                               Field("\\chi" + str(type),
                                     [], [2], dim + half),
                               2 * Field("c" + str(type),
                                     dim = dim + 1) ],
                              [Field("\\psi" + str(type),
                                     [1], [], dim + half),
                               Field("v" + str(type),
                                     [1], [1], dim + 1),
                               Field("v" + str(type),
                                     [1], [2], dim + 1),
                               2 * Field("\\eta" + str(type),
                                     [1], [], dim + 1 + half) ],
                              [Field("\\psi" + str(type),
                                     [2], [], dim + half),
                               Field("v" + str(type),
                                     [2], [1], dim + 1),
                               Field("v" + str(type),
                                     [2], [2], dim + 1),
                               2 * Field("\\eta" + str(type),
                                     [2], [], dim + 1 + half) ],
                              [-2 * Field("b" + str(type),
                                     dim = dim + 1),
                               -2 * Field("\\xi" + str(type),
                                     [], [1], dim + 1 + half),
                               -2 * Field("\\xi" + str(type),
                                     [], [2], dim + 1 + half),
                               -4 * Field("d" + str(type),
                                     dim = dim + 2)]]

    def copy(self):
        """Creates a copy of the supermultiplet
        """
        other = SuperField(dim = self.dim,
                           Grassmann = self.Grassmann)
        other.component = [[field.copy()
                            for field in row]
                           for row in self.component]
        return other

    def subs(self, x, y):
        """Returns a supermultiplet with x substituted by y
        """
        other = SuperField(dim = self.dim,
                           Grassmann = self.Grassmann)
        other.component = [[field.subs(x, y)
                            for field in row]
                           for row in self.component]
        return other

    def check(self):
        """Additional checks about the scaling dimensions"""
        # (1) check that the superfield is indeed a 4x4 matrix
        if len(self.component) != 4:
            raise Error("Superfield::check()",
                        "supermultiplet does noe have 4 rows",
                        "")
        for row in self.component:
            if len(row) != 4:
                raise Error("Superfield::check()",
                            "supermultiplet does noe have 4 columns",
                            "")
        # (2) check all components
        for row in self.component:
            for field in row:
                field.check()
        # (3) check that the scaling dimensions match
        dim = self.component[0][0].dim

        t = tuple([ tuple([field.dim - dim for field in row])
                    for row in self.component])
        if t != dimTable:
            raise Error("Superfield::check()",
                        "incorrect scaling dimensions",
                        "")
        # (4) check that the spin statistics match
        Grassmann = self.component[0][0].Grassmann
        t = tuple([ tuple([field.Grassmann != Grassmann
                           for field in row])
                    for row in self.component])
        if t != GrassmannTable:
            raise Error("Superfield::check()",
                        "incorrect spin statistics",
                        "")

    def printlowest(self, latexoutput=latexoutput):
        """Print the lowest element of the multiplet
        """
        self.component[0][0].print(latexoutput)

    def printall(self):
        """Print all 16 components of the multiplet
        """
        for i in range(4):
            for j in range(4):
                self.component[i][j].print()

    def lowest(self):
        """Returns the lowest component of the superfield
        """
        return self.component[0][0]

    def Fterm(self):
        """Returns the F term, that is the theta^2 component of the superfield
        """
        return -half * self.component[3][0]

    def Fbarterm(self):
        """Returns the F^bar term, that is the (theta^bar)^2 component of the superfield
        """
        return half * self.component[0][3]

    def Dterm(self):
        """Return the D term, that is the (theta theta^bar)^2 component
        """
        return -quarter * self.component[3][3]

    def iszero(self):
        """Check if a superfield is zero
        """
        for i in range(4):
            for j in range(4):
                if self.component[i][j].iszero() == False:
                    return False
        return True

    def sort(self):
        """Simplifies each component of the superfield
        """
        for row in self.component:
            for field in row:
                field.sort()

    def d(self, mu):
        """derivative acting on the expression
        """
        other = SuperField(dim = self.dim + 1,
                           Grassmann=self.Grassmann)
        other.component = [ [field.d(mu)
                             for field in row]
                            for row in self.component]
        return other
        
    def dd(self, a, b):
        """Alternate definition of the derivative using spin indices
        """
        return sum(
            lambda aa, bb, mu:
                sigma[mu[0]][a][b],
            lambda aa, bb, mu:
                self.d(mu[0]),
            0, 0, 1
        )

    def box(self):
        """Laplacian of the superfield
        """
        return (self.d(0).d(0) - self.d(1).d(1)
                - self.d(2).d(2) - self.d(3).d(3))

    def __add__(self, other):
        """Sum of two supermultiplets
        """
        if self.dim != other.dim:
            raise Error("SuperField::__add__",
                        "Adding superfields with different "
                        "scaling dimensions", "")
        if self.Grassmann != other.Grassmann:
            raise Error("SuperField::__add__",
                        "Adding superfields with different "
                        "spin statistics", "")
        result = SuperField(dim = self.dim,
                            Grassmann=self.Grassmann)
        result.component = [ [self.component[i][j]
                              + other.component[i][j]
                              for j in range(4)]
                             for i in range(4)]
        return result

    def __rmul__(self, c):
        """Multiplication by a real number
        """
        result = SuperField(dim = self.dim,
                            Grassmann=self.Grassmann)
        result.component = [ [c * field
                              for field in row]
                             for row in self.component]
        return result

    def __sub__(self, other):
        """Difference of two multiplets
        """
        return self + (-1) * other

    def __mul__(self, other):
        """Multiplication of two superfields
        """
        result = SuperField(dim = self.dim + other.dim,
                            Grassmann = (self.Grassmann != other.Grassmann))
        if other.Grassmann:
            Rule = MultGrassmannRule
        else:
            Rule = MultRule
        for rule in Rule:
            result.component[rule[1]][rule[2]] += (
                rule[0]
                * self.component[rule[3]][rule[4]]
                * other.component[rule[5]][rule[6]]
            )
        return result

    def covDcomputation(self, Rule):
        """Computation of the covariant derivatives
        """
        result = SuperField(dim = self.dim + half,
                            Grassmann = not self.Grassmann)
        table = ([self.component]
             + [ [ [ field.d(mu) for field in row ]
                   for row in self.component ]
                 for mu in range(4)])
        for rule in Rule:
            result.component[rule[1]][rule[2]] += (
                rule[0] * table[rule[3]][rule[4]][rule[5]]
            )
        result.sort()
        if self.Grassmann:
            result = (-1) * result
        return result

    def covD(self, a):
        """Covariant derivative with respect to theta
        """
        return self.covDcomputation(covDRule[a])
    def covDbar(self, a):
        """Covariant derivative with respect to theta^bar
        """
        return self.covDcomputation(covDbarRule[a])

    def covD2(self):
        return 2 * self.covD(0).covD(1)
    def covDbar2(self):
        return 2 * self.covDbar(1).covDbar(0)

    def coefficientset(self):
        """Set of coefficients of the terms"""
        coeffs = set()
        for row in self.component:
            for field in row:
                coeffs.update(field.coefficientset())
        return coeffs

    def linearcombinationcoeffs(self, superfieldlist, variables = []):
        """Returns the coefficients of the linear combination
        of a list of superfields giving the present superfield
        (note 1: the result can be [] )
        (note 2: one can specify additional variables)
        """
        n = len(superfieldlist)
        if n == 0:
            return []
        # (1) construct a linear combination X of all the superfields
        #     with symbolic coefficients p
        p0 = []
        X = self.copy()
        for superfield in superfieldlist:
            q = symbols("ppp" + str(len(p0)))
            p0.append(q)
            X -= q * superfield
        X.sort()
        p = p0 + variables
        # (2) try to find a solution to the equation X = 0
        # by writing it in the form A * p = B
        Xc = list(X.coefficientset())
        A = [[diff(XX, pp) for pp in p]
             for XX in Xc]
        B = [XX.subs([(pp,0) for pp in p]) for XX in Xc]
        sol = linsolve( (Matrix(A), -Matrix(B)), p)
        # check the existence of solutions
        if len(sol) == 0:
            if len(variables) > 0:
                raise Error("Superfield::linearcombinationcoeffs",
                            "no solution", "")
            return []
        if len(sol) > 1:
            raise Error("Superfield::linearcombinationcoeffs",
                        "multiple solutions", "")
        sol = list(sol)[0]
        # check that the result can be written in terms of the variables only
        if len(set(p0).intersection(sol.free_symbols)) > 0:
            raise Error("Superfield::linearcombinationcoeffs",
                        "degenerate set of fields in argument", "")
        # if there are additional variables, return the solution for them only
        if len(variables) > 0:
            return [p[i] - sol[i] for i in range(n,len(p))
                    if p[i] not in sol.free_symbols]
        # if there are no additional variables, return the coefficients
        return sol

    def islinearcombination(self, superfieldset):
        """Check if the superfield is a linear combination of a set of terms
        """
        x = self.linearcombinationcoeffs(superfieldset)
        if len(x) > 0:
            return True
        return False

def reducebasis(basis):
    """Reduces a basis to only linearly independent terms
    """
    n = len(basis)
    print("\nReducing a basis of", n, "elements:", flush=True)
    result = []
    i = 0
    for element in basis:
        i += 1
        print(" -> analyzing " + str(i) + "/" + str(n),
              end="", flush=True)
        element.sort()
        if not element.islinearcombination(result):
            result.append(element)
            print("", flush=True)
        else:
            print("  -> NOT linearly independent!", flush=True)
    print(" ->", len(result),
          "are linearly independent.\n", flush=True)
    return result

##############################################################################


def printbegin(expression):
    """Additional function used for displaying results
    """
    print("[{:%H:%M:%S}] >>>".format(datetime.datetime.now()),
          expression, end="", flush=True)
def printend(expression = "ok"):
    """Additional function used for displaying results
    """
    print(expression, flush=True)
