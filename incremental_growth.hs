import Helper

type Alpha = (Double, Double) -> Double -> ((Double, Double) -> Double -> Double) -> Double -> Double
type Growth  = (Double, Double) -> Double -> Double

data State = State 
            { bounds :: (Double, Double)
            , gCum :: Growth
            , alphaCum :: Alpha
            }

instance Show State where
    show state = "State { bounds = " ++ show (bounds state) ++" }"

initialState = State {
        bounds = (1,2),
        gCum = g,
        alphaCum = alpha
    }    

mu = 100
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

r :: (Double, Double) -- lower bound in reference configuration
  -> Double -- point in reference configuration
  -> ((Double, Double) -> Double -> Double) -- g
  -> Double -- a
  -> Double
r bounds x g a = 
    let integrand s = s**2 * (g bounds s)**3
        lowerBound = fst bounds
    in (a**3 + 3*(quadrature64 integrand lowerBound x))**(1/3)

alpha :: (Double, Double) --lowerBound
      -> Double -- x
      -> ((Double, Double) -> Double -> Double) -- g
      -> Double -- a
      -> Double
alpha bounds x g a = (r bounds x g a)/(x * (g bounds x))

t :: (Double, Double) -- lower bound in reference configuration
  -> Double -- point in reference configuration
  -> ((Double, Double) -> Double -> Double) -- g
  -> ((Double, Double) -> Double -> ((Double, Double)  -> Double -> Double) -> Double -> Double) -- alpha
  -> Double -- lower bound in current configuration
  -> Double -- stress at a point in current configuration
t bounds x g alpha a = 
    let integrand' s = 2*mu/s*((alpha bounds s g a)**(-1) - (alpha bounds s g a)**(-7))
        lowerBound = fst bounds
    in quadrature64 integrand' lowerBound x

updateState :: State -> State
updateState state = 
    let upperBound = snd (bounds state)
        alpha = alphaCum state
        g = gCum state
        a = secant 1e-14 (t (bounds state) upperBound g alpha) (bounds state)
    in state {
        bounds = (a, r (bounds state) upperBound g a),
        gCum = \s x -> (g (bounds state) upperBound) * g s x,
        alphaCum = \x1 x2 x3 x4 -> (alpha (bounds state) upperBound g a) * ((alpha x1 x2 x3 x4))
    }

fst' :: (a, b, c) -> a
fst' (a, _, _) = a

volume :: (Double, Double) -> Double 
volume (inner, outer) = 4/3*pi*(outer^3 - inner^3)

volumeRatio :: (Double, Double) -> (Double, Double) -> Double
volumeRatio (inner, outer) (innerNew, outerNew) = (volume (innerNew, outerNew)) / (volume (inner, outer))

main :: IO ()
main = do
    let initialState = State {
        bounds = (1,2),
        gCum = g,
        alphaCum = alpha
    }    
    let tenSteps = take 9 $ iterate updateState initialState

    print $ map bounds tenSteps