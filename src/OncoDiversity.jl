module OncoDiversity

using DataFrames
using CSV
using Query
using Statistics
using Dates
using Plots
import DataFrames: rename!

export diversity, calcOverQ, findInflection, findInflectionLocal, reportInflection, calcIPdiversity, patientdiversity,
  rename!, requirecol, countframe, countframefull, summarize, filteredDataset, 
  Dataset, CDR3SeqData,
  AbstractAnalysis, CDR3SeqAnalysis
  

abstract type AbstractAnalysis end
struct CDR3SeqAnalysis <: AbstractAnalysis end
struct SCRNASeqAnalysis<: AbstractAnalysis end
struct FlowCytoAnalysis <: AbstractAnalysis end
struct CytokineAnalysis <: AbstractAnalysis end
struct SomaScanAnalysis <: AbstractAnalysis end

abstract type Dataset end

struct CDR3SeqData{T,U} <: Dataset
  counts::T
  physiochem::U
end

function diversity(seq, q)
  0 <= q || error("q=$(q) must be greater than zero")

  if(q == 1.0)
    return exp(-sum(seq .* log.(seq)))
  else
    return sum(seq.^q)^(1/(1-q))
  end

end 

function calcOverQ(origDF::DataFrame,rangeQ::Array{Float64})
  qDdf = DataFrame(patient=unique(origDF[:,:patient]))

  for q in rangeQ
      tempDF = patientdiversity(origDF,q);
      
      qDdf[:,Symbol(q)]=tempDF[:,:diversity]
  end
  
  return qDdf
end

function findInflection(q::Array, qD::Array)
  approxLogDeriv = q[2:end].*diff(qD)./diff(q);
  return approxLogDeriv; 
end

function findInflectionLocal(approx::Array)
  return argmin((approx))
end

function reportInflection(q::Array, qD::Array)
  diffVector = findInflection(q, qD)
  inflectPt = findInflectionLocal(diffVector)
  slopeInflectPt = abs(diffVector[inflectPt])
  qInflecPt = q[inflectPt+1]
  return (qInflecPt, slopeInflectPt)
end

function calcIPdiversity(df::DataFrame, qrange::Array)
    
  # set up new dataframe to capture the inflection point (InfPt) and slope at the inflection point (InflPtS)
  newDF = DataFrame(
            patient=unique(df[:,:patient]),
            InflPt=zeros(size(unique(df[:,:patient]))[1]),
            InflPtS=zeros(size(unique(df[:,:patient]))[1]))
  
  # calculate qD over a supplied range of q's
  dfRangeQ = calcOverQ(df,qrange)
  
  # calculate the IP & slope and fill into the dataframe
  counter = 1; 
  for r in eachrow(dfRangeQ)
       
      qDs = collect(r[2:end])

      (IP,IPslope) = reportInflection(qrange,qDs)
      newDF[counter,:InflPt] = IP
      newDF[counter,:InflPtS] = IPslope
      counter = counter+1; 
  end
  
  return newDF
end

function patientdiversity(countframe::DataFrame,q::Real)
  #@show totUniq = length(unique(df.CDR3))
      
  diversityframe = @from i in countframe begin
          @group i by i.patient into g
          @select {
              patient=first(g.patient), 
              diversity=diversity(g.count ./ sum(g.count), q),
              numclust = length(g.count),
              checkrange = length(g.count) > diversity(g.count ./ sum(g.count), q)
          }
      @collect DataFrame
  end
  return diversityframe
end

function requirecol(df::DataFrame, colname::AbstractString) 
  @assert colname in names(df) "Expected column named: $colname, names are $(names(df))"
end

rename!(df::DataFrame, s::String, t::Symbol) = rename!(df, Symbol(s)=>t)

function countframe(a::CDR3SeqAnalysis, typeDF::DataFrame)
  requirecol(typeDF, "patient")
  requirecol(typeDF, "CDR3")
  
  countDF = @from i in typeDF begin
    @group i by i.patient, i.CDR3 into g
    @select {
        patient=first(g.patient), 
        CDR3=first(g.CDR3), 
        count=length(g.CDR3)
    }
    @collect DataFrame
  end;
end

function countframefull(a::CDR3SeqAnalysis, df::DataFrame)
  rename!(df, "Filename", :patient)
  return innerjoin(countframe(CDR3SeqAnalysis(), df), df, on=[:patient, :CDR3])
end


function CDR3SeqData(df::DataFrame)
  requirecol(df, "CDR3")
  try
    requirecol(df, "patient")
  catch AssertionError
    requirecol(df, "Filename")
    rename!(df, "Filename", :patient)
  end
  return CDR3SeqData(
    countframe(CDR3SeqAnalysis(), df),
    select(df, Not(:patient)) |> @groupby(_.CDR3) |> x-> map(first, x) |> DataFrame
  )
end

DataFrame(data::CDR3SeqData) = innerjoin(data.counts, data.physiochem, on=[:CDR3])

function summarize(df::CDR3SeqData)
  statframe = @from i in df.counts begin
  @group i by i.patient into g
  @select {
      patient=first(g.patient), 
      total_count=sum(g.count), 
      min_count = minimum(g.count),
      max_count = maximum(g.count),
      avg_count = mean(g.count),
      median_count = median(g.count),
      num_unique_CDR3 = length(g.CDR3)
  }
  @collect DataFrame
  end
      
  return statframe
      
end


function filteredDataset(::Type{CDR3SeqData}, df::DataFrame, receptorname::String)
  ## subset the type group only 
  recepframe = @from i in df begin
    @where i.Receptor == receptorname
    @select i
    @collect DataFrame
  end;
  println("Number of recoveries: ", size(recepframe))
  println("Cohort size: ", length(unique(recepframe[:,:Filename])))
  ds = CDR3SeqData(recepframe)
  return ds
end

# function foo(c::CDR3SeqAnalysis, df::DataFrame)
#   #...CDR3Analysis code
# end

# function foo(c::SCRNASeqAnalysis, df::DataFrame)
#   #...SCRNASeq code here
# end

end
