import Helper 

data State = State 
            { innerBoundary :: Double   -- inner boundary
            , outerBoundary :: Double   -- outer boundary
            , gammaR :: Double          -- radial growth rate
            , gammaTheta :: Double      -- angular growth rate
            , alphaR :: Double          -- radial strain coefficient 
            , alphaTheta :: Double      -- angular strain coefficient
            , mu :: Double              -- stress/strain ratio (Young modulus)
            } deriving (Show)

-- r = (a^2 + γr γθ (R2 − A2))^(1/2)
updateR :: State
        -> Double   -- inner radius in new configuration          
        -> Double   -- point in reference (previous) configuration
        -> Double   -- r
updateR st a oldR = sqrt(a ** 2 + (gammaR st) * (gammaTheta st) * (oldR ** 2 - (innerBoundary st) ** 2))

-- αr^(i+1) = αr^i R γθ / r
updateAlphaR :: State 
             -> Double  -- inner radius in new configuration          
             -> Double  -- point in reference (previous) configuration
             -> Double  -- αr^(i+1)
updateAlphaR st a oldR = 
    let r = updateR st a oldR
    in  (alphaR st) * oldR * (gammaTheta st) / r

-- αθ^(i+1) = αθ^i r / (R γθ)
updateAlphaTheta :: State 
                 -> Double  -- inner radius in new configuration          
                 -> Double  -- point in reference (previous) configuration
                 -> Double  -- αθ^(i+1)
updateAlphaTheta st a oldR = 
    let r = updateR st a oldR
    in  (alphaTheta st) * r / (oldR * (gammaTheta st))

-- P = μγr γθ ∫_A^B s r^2 (α_θ)^2 − (α_r)^2 ds
calculateP :: State 
           -> Double    -- inner radius in new configuration
           -> Double    -- P
calculateP st a = 
    let r = updateR st a 
        alphaR_ = updateAlphaR st a  
        alphaTheta_ = updateAlphaTheta st a 

        integrand = (\s -> s/(r s)**2 * ((alphaTheta_ s) ** 2 - (alphaR_ s) ** 2))
        integral = quadrature64 integrand (innerBoundary st) (outerBoundary st)

    in (mu st) * (gammaR st) * (gammaTheta st) * integral

findInnerBoundary :: State -> Double
findInnerBoundary st = secant 1e-14 (\a -> calculateP st a) ((innerBoundary st), (outerBoundary st))  -- here i am assuming the new inner boundary is between the old boundaries

-- σrr = μγr γθ ∫_A^R s r^2 (α_θ)^2 − (α_r)^2 ds − P
sigmaRR :: State 
        -> Double   -- inner radius in new configuration          
        -> Double   -- point in reference (previous) configuration
        -> Double   -- σrr
sigmaRR st a oldR = 
    let r = updateR st a 
        alphaR_ = updateAlphaR st a  
        alphaTheta_ = updateAlphaTheta st a 
        p = calculateP st a

        integrand = (\s -> s / (r s)**2 * ((alphaTheta_ s) ** 2 - (alphaR_ s) ** 2))
        integral = quadrature64 integrand (innerBoundary st) oldR
    in  (mu st) * (gammaR st) * (gammaTheta st) * integral - p

-- σθθ = μγr γθ ∫_A^R s r^2 (α_θ)^2 − (α_r)^2 ds + σrr
sigmaTT :: State 
        -> Double   -- inner radius in new configuration          
        -> Double   -- point in reference (previous) configuration 
        -> Double   -- σθθ
sigmaTT st a oldR = 
    let r = updateR st a 
        alphaR_ = updateAlphaR st a  
        alphaTheta_ = updateAlphaTheta st a 

        dSigmaRR = (\s -> s / (r s)**2 * ((alphaTheta_ s) ** 2 - (alphaR_ s) ** 2))
    in  (mu st) * (gammaR st) * (gammaTheta st) * (dSigmaRR oldR) * (r oldR)**2 / (oldR * (gammaR st) * (gammaTheta st)) + (sigmaRR st a oldR)

main :: IO ()
main = do
    let initialState = State {
        innerBoundary = 1.0,
        outerBoundary = 2.0,
        gammaR = 1,
        gammaTheta = 1.5,
        alphaR = 1,
        alphaTheta = 1,
        mu = 1
    }    
    let a = findInnerBoundary initialState 
        radialStress = [sigmaRR initialState a r | r <- equidistant 10 1 2]
    print $ radialStress