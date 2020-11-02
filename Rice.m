classdef Rice
    %RICE This class holds all the parameters for the Rician fading model.
    % It calculates theoretical envelope and phase PDFs
    % It does a Monte Carlo simulation using the parameters
    
    properties(Constant, Hidden = true)
        NumSamples = 2E6; % number of samples
        r = 0:0.001:6 % envelope range for PDF ploteter
        theta = -pi:0.001:pi; % phase range for PDF plot
    end 
    
    properties(Access = public)
        K; % Rician K factor
        r_hat_2; % root mean square of the signal
        phi; % phase parameter
    end
    
    properties(Hidden = true) 
        multipathFading; % Found based on the inputs
        envelopeProbability; % Calculated theoretical envelope probability
        phaseProbability; % Calculated theoretical phase probability
        xdataEnv; % Simulated envelope density plot x values 
        ydataEnv; % Simlated envelope density plot y valyes
        xdataPh; % Simulated phase density plot x values 
        ydataPh; % Simlated phase density plot y valyes
    end
    
    methods(Access = public)
        function obj = Rice(K,r_hat_2,phi)
            %ADDITIVESHADOWRICE Construct an instance of this class
            %   Assigning input values
            
            % K must be positive valued and restricting to a maximum of 50
            if ~isnumeric(K)
               error('Error. \nK must be a number, not a %s.', class(K));
            end
            
            if K < 0 || K > 50
               error('Error. \n%s','K must be greater than 0 and less than 50.');
            end
         
            obj.K = K;
            
            % r_hat_2 must be greater than 0.5 and less than 2.5
            if ~isnumeric(r_hat_2)
               error('Error. \nr_hat_2 must be a number, not a %s.', class(K));
            end
            
            if r_hat_2 < 0.5 || r_hat_2 > 2.5
               error('Error. \n%s','r_hat_2 must be greater than 0.5 and less than 2.5.');
            end
            
            obj.r_hat_2 = r_hat_2;
            
            % phi is in the range -pi to pi
            if ~isnumeric(phi)
               error('Error. \nphi must be a number, not a %s.', class(K));
            end
            
            if phi < -pi || phi > pi
               error('Error. \n%s','phi must be between -pi and pi.');
            end
            
            obj.phi = phi;
            
            % other calculated properties
            obj.multipathFading = complexMultipathFading(obj);
            obj.envelopeProbability = envelopePDF(obj);
            obj.phaseProbability = phasePDF(obj);
            [obj.xdataEnv, obj.ydataEnv] = envelopeDensity(obj);
            [obj.xdataPh, obj.ydataPh] = phaseDensity(obj);
            
        end
    end
    
    methods(Access = private)
        function [p, q] = means(obj)
            %means Calculates the means of the complex Gaussians 
            %representing the in-phase and quadrature components.

            p = sqrt(obj.K/(1+obj.K)) .* sqrt(obj.r_hat_2) .* cos(obj.phi);
            q = sqrt(obj.K/(1+obj.K)) .* sqrt(obj.r_hat_2) .* sin(obj.phi);

        end
        
        function [sigma] = scatteredComponent(obj)
            %scatteredComponent Calculates the power of the scattered 
            %signal component.    
            
            sigma = sqrt(obj.r_hat_2 ./(2 * (1 + obj.K)));
        
        end
        
        function [gaussians] = generateGaussians(obj, mean, sigma) 
            %generateGaussians Generates the Gaussian random variables 
            
            gaussians = normrnd(mean,sigma,[1,obj.NumSamples]);
        end
        
        function [multipathFading] = complexMultipathFading(obj) 
            %complexMultipathFading Generates the Rician fading model 
            
            [p, q] = means(obj);
            [sigma] = scatteredComponent(obj);
            
            multipathFading = generateGaussians(obj, p, sigma) + 1i.* generateGaussians(obj, q, sigma);
        end    
        
        function [eProbTheor] = envelopePDF(obj)
            %envelopePDF Calculates the theoretical envelope PDF
            
            eProbTheor = 2 .* (1+obj.K) .* obj.r ./ obj.r_hat_2 ...
                .* exp(-obj.K - ((1+obj.K).*obj.r.^2 ./ obj.r_hat_2))...
                .* besseli(0, 2.*obj.r.*sqrt(obj.K.*(obj.K+1)./obj.r_hat_2));
        end
        
        function [pProbTheor] = phasePDF(obj)
            %envelopePDF Calculates the theoretical phase PDF
            
            pProbTheor = 1./(2*pi).*exp(-obj.K) .*(1 + (sqrt(4.*pi.*obj.K)...
                .* exp(obj.K.*(cos(obj.theta-obj.phi)).^2) .* cos(obj.theta-obj.phi))...
                .* (1 - qfunc(sqrt(2.*obj.K).*cos(obj.theta-obj.phi))));
        end
        
        function [xdataEnv, ydataEnv] = envelopeDensity(obj)
            %envelopeDensity Evaluates the envelope PDF
            R = sqrt((real(obj.multipathFading)).^(2) + (imag(obj.multipathFading)).^(2));

            [f,x] = ecdf(R);
            [ydataEnv, xdataEnv] = ecdfhist(f,x, 0:0.05:max(obj.r));
        end
        
        function [xdataPh, ydataPh] = phaseDensity(obj)
            %envelopeDensity Evaluates the envelope PDF
            T = angle((real(obj.multipathFading)) + (imag(obj.multipathFading)));

            [f,x] = ecdf(T);
            [ydataPh, xdataPh] = ecdfhist(f,x, min(obj.theta)+0.03:0.06:max(obj.theta));
        end
            
    end
end

