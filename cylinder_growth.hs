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
updateR state a oldR =
    let gammaR' = gammaR state
        gammaTheta' = gammaTheta state
        innerBoundary' = innerBoundary state
    in  sqrt(a ** 2 + gammaR' * gammaTheta' * (oldR ** 2 - innerBoundary' ** 2))

-- αr^(i+1) = αr^i R γθ / r
updateAlphaR :: State 
             -> Double  -- inner radius in new configuration          
             -> Double  -- point in reference (previous) configuration
             -> Double  -- αr^(i+1)
updateAlphaR state a oldR = 
    let gammaTheta' = gammaTheta state
        alphaR' = alphaR state
        r' = updateR state a oldR
    in  alphaR' * oldR * gammaTheta' / r'

-- αθ^(i+1) = αθ^i r / (R γθ)
updateAlphaTheta :: State 
                 -> Double  -- inner radius in new configuration          
                 -> Double  -- point in reference (previous) configuration
                 -> Double  -- αθ^(i+1)
updateAlphaTheta state a oldR = 
    let gammaTheta' = gammaTheta state
        alphaTheta' = alphaTheta state
        r' = updateR state a oldR
    in  alphaTheta' * r' / (oldR * gammaTheta')

-- P = μγr γθ ∫_A^B s r^2 (α_θ)^2 − (α_r)^2 ds
calculateP :: State 
           -> Double    -- inner radius in new configuration
           -> Double    -- P
calculateP state a = 
    let mu' = mu state
        gammaR' = gammaR state
        gammaTheta' = gammaTheta state
        alphaTheta' = alphaTheta state
        alphaR' = alphaR state
        innerBoundary' = innerBoundary state
        outerBoundary' = outerBoundary state

        r = (\s -> sqrt(a ** 2 + gammaR' * gammaTheta' * (s ** 2 - innerBoundary' ** 2)))
        alphaR_ = (\s -> alphaR' * s * gammaTheta' / (r s))
        alphaTheta_ = (\s -> alphaTheta' * (r s) / (s * gammaTheta'))

        integrand = (\s -> s/(r s)**2 * ((alphaTheta_ s) ** 2 - (alphaR_ s) ** 2))
        integral = quadrature64 integrand innerBoundary' outerBoundary'
    in mu' * gammaR' * gammaTheta' * integral

findInnerBoundary :: State -> Double
findInnerBoundary state = 
    let innerBoundary' = innerBoundary state
        outerBoundary' = outerBoundary state
    in  secant 1e-14 (\a -> calculateP state a) (innerBoundary', outerBoundary')  -- here i am assuming the new inner boundary is between the old boundaries

-- σrr = μγr γθ ∫_A^R s r^2 (α_θ)^2 − (α_r)^2 ds − P
sigmaRR :: State 
        -> Double   -- inner radius in new configuration          
        -> Double   -- point in reference (previous) configuration
        -> Double   -- σrr
sigmaRR state a oldR = 
    let mu' = mu state
        gammaR' = gammaR state
        gammaTheta' = gammaTheta state
        alphaTheta' = alphaTheta state
        alphaR' = alphaR state
        innerBoundary' = innerBoundary state

        r = (\s -> sqrt(a ** 2 + gammaR' * gammaTheta' * (s ** 2 - innerBoundary' ** 2)))
        alphaR_ = (\s -> alphaR' * s * gammaTheta' / (r s))
        alphaTheta_ = (\s -> alphaTheta' * (r s) / (s * gammaTheta'))
        p = calculateP state a

        integrand = (\s -> s / (r s)**2 * ((alphaTheta_ s) ** 2 - (alphaR_ s) ** 2))
        integral = quadrature64 integrand innerBoundary' oldR
    in  mu' * gammaR' * gammaTheta' * integral - p

-- σθθ = μγr γθ ∫_A^R s r^2 (α_θ)^2 − (α_r)^2 ds + σrr
sigmaTT :: State 
        -> Double   -- inner radius in new configuration          
        -> Double   -- point in reference (previous) configuration 
        -> Double   -- σθθ
sigmaTT state a oldR = 
    let mu' = mu state
        gammaR' = gammaR state
        gammaTheta' = gammaTheta state
        alphaTheta' = alphaTheta state
        alphaR' = alphaR state
        innerBoundary' = innerBoundary state
        outerBoundary' = outerBoundary state

        r = (\s -> sqrt(a ** 2 + gammaR' * gammaTheta' * (s ** 2 - innerBoundary' ** 2)))
        alphaR_ = (\s -> alphaR' * s * gammaTheta' / (r s))
        alphaTheta_ = (\s -> alphaTheta' * (r s) / (s * gammaTheta'))

        dSigmaRR = (\s -> s / (r s)**2 * ((alphaTheta_ s) ** 2 - (alphaR_ s) ** 2))
    in  mu' * gammaR' * gammaTheta' * (dSigmaRR oldR) * (r oldR)**2 / (oldR * gammaR' * gammaTheta') + (sigmaRR state a oldR)

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