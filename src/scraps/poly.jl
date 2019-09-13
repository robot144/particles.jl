#struct Poly{T}
#  a::Vector{T}
#end


function add{T}(b1::Poly{T},b2::Poly{T})
  if length(b1.a)>=length(b2.a)
      temp=copy(b1.a)
      for i=1:length(b2.a)
         temp[i]+=b2.a[i]
      end
   else
      temp=copy(b2.a)
      for i=1:length(b1.a)
         temp[i]+=b1.a[i]
      end
   end
   return Poly(temp)
end

function display{T}(a::Poly{T})
   println("TODO")
end

a1=Poly([1.0,2.0])
a2=Poly([1.1,2.1,3.1])

b=add(a1,a2)
