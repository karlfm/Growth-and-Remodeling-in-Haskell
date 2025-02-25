import Helper

type Alpha = (Double, Double) -> Double -> ((Double, Double) -> Double -> Double) -> Double -> Double
type Growth  = (Double, Double) -> Double -> Double

data State = State 
            { bounds :: (Double, Double)
            , gCum :: Growth
            , alphaCum :: Alpha
            }

instance Show State where
    show state = "State { bounds = " ++ show (bounds state) ++ " }"

mu = 100 -- Shear modulus

growthRate :: (Double, Double) 
           -> Double
growthRate bounds = 
    let a = fst bounds
        b = snd bounds
    in (((b**3 - a**3) * 1.01 / b**3 + (a/b)**3)**(1/3) - 1)/(b - a)

g :: (Double, Double) -- inner radius in reference configuration
  -> Double -- point in reference configuration
  -> Double -- point in current configuration
g bounds x = 1 + (growthRate bounds) * (x - (fst bounds))

-- radius in spatial configuration
r :: (Double, Double) -- lower bound in reference configuration
  -> Double -- point in reference configuration
  -> Growth
  -> Double -- a
  -> Double
r bounds x g a = 
    let integrand s = s**2 * (g bounds s)**3
        lowerBound = fst bounds
    in (a**3 + 3*(quadrature64 integrand lowerBound x))**(1/3)

-- elastic strain
alpha :: (Double, Double) --lowerBound
      -> Double -- x
      -> Growth
      -> Double -- a
      -> Double
alpha bounds x g a = (r bounds x g a)/(x * (g bounds x))

-- stress in spatial configuration (Cauchy stress)
t :: (Double, Double) -- lower bound in reference configuration
  -> Double -- point in reference configuration
  -> Growth
  -> Alpha
  -> Double -- lower bound in current configuration
  -> Double -- stress at a point in current configuration
t bounds x g alpha a = 
    let integrand s = 2*mu/s*((alpha bounds s g a)**(-1) - (alpha bounds s g a)**(-7))
        lowerBound = fst bounds
    in quadrature64 integrand lowerBound x

updateState :: State -> State
updateState state = 
    let upperBound = snd (bounds state)
        alpha = alphaCum state
        g = gCum state
        a = secant 1e-14 (t (bounds state) upperBound g alpha) (bounds state)
        oldGCum = g (bounds state) upperBound
        oldAlphaCum = alpha (bounds state) upperBound g a
    in state {
        bounds = (a, r (bounds state) upperBound g a),
        gCum = \s x -> oldGCum * g s x,
        alphaCum = \x1 x2 x3 x4 -> oldAlphaCum * (alpha x1 x2 x3 x4)
    }

main :: IO ()
main = do
    let initialState = State {
        bounds = (1,2),
        gCum = g,
        alphaCum = alpha
    }    
    let steps = take 6 $ iterate updateState initialState

    print $ map bounds steps