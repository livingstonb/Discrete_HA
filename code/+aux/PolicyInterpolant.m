classdef PolicyInterpolant < handle
	properties
		dims;
		n_dims;

		gridded_interpolant;
		xmins;
		ymins;
	end

	methods
		function obj = PolicyInterpolant(X, Y, varargin)
			shape = size(X);
			tmp = shape(2:end);
			
			obj.n_dims = numel(tmp(:));
			obj.dims = [tmp(:); ones(5-obj.n_dims, 1)];
			obj.xmins = squeeze(X(1,:,:,:,:,:));
			obj.ymins = squeeze(Y(1,:,:,:,:,:));
		end

		function out = interp(obj, qvals, varargin)
			out = obj.gridded_interpolant{varargin{:}}(qvals);
		end

		function out = interp_extended_below(obj, qvals, varargin)
			out = size(qvals);

			adjust = qvals < obj.xmins(varargin{:});
			out(~adjust) = obj.gridded_interpolant{varargin{:}}(qvals(~adjust));
			out(adjust) = obj.ymins(varargin{:});
		end

		function set_interpolant(obj, x, y, varargin)
			obj.gridded_interpolant{varargin{:}} = ...
				griddedInterpolant(x(:), y(:));
		end
	end

	methods (Access=private)
		

		function loop_over_dims(obj, fn_handle, fn_args)
			for i5 = 1:obj.dims(5)
				for i4 = 1:obj.dims(4)
					for i3 = 1:obj.dims(3)
						for i2 = 1:obj.dims(2)
							for i1 = 1:obj.dims(1)
								idims = {i1, i2, i3, i4, i5};
								fn_handle(idims{:}, fn_args{:});
							end
						end
					end
				end
			end
		end
	end
end