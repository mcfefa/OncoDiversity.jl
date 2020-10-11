module OncoDiversity

using DataFrames
using CSV
using Query
using Statistics
using Dates
import DataFrames: rename!

export diversity, calcOverQ, findInflection, findInflectionLocal, reportInflection, calcIPdiversity, patientdiversity,
  rename!, requirecol, countframe, countframefull, summarize, filteredDataset, 
  Dataset, CDR3SeqData,
  AbstractAnalysis, CDR3SeqAnalysis,
  DiversityScores
  

abstract type AbstractAnalysis end
struct CDR3SeqAnalysis <: AbstractAnalysis end
struct SCRNASeqAnalysis<: AbstractAnalysis end
struct FlowCytoAnalysis <: AbstractAnalysis end
struct CytokineAnalysis <: AbstractAnalysis end
struct SomaScanAnalysis <: AbstractAnalysis end

abstract type Dataset end

"""    CDR3SeqData{T,U}

store two data tables that make up a CDR3 Dataset,
the per (patient, sequence) count and the physiochemisty data for each sequence.
These tables can be joined together usind `DataFrame(ds::CDR3SeqData)`.

"""
struct CDR3SeqData{T,U} <: Dataset
  counts::T
  physiochem::U
end

"""    diversity(seq, q)

Computes the generalized diversity index for a given value of q and the associated frequencies of the dataset. 

```math
  ^qD = \\left( \\sum_i^n (p_i^q) \\right)^{1/(1-q)}
```
where p_i are the frequency across types and q is the order of diversity and you are summing over all types present in the dataset.

"""
function diversity(seq, q)
  0 <= q || error("q=$(q) must be greater than zero")

  if(q == 1.0)
    return exp(-sum(seq .* log.(seq)))
  else
    return sum(seq.^q)^(1/(1-q))
  end

end 

"""    calcOverQ(origDF::DataFrame,rangeQ::Array{Float64})

Computes the diversity score over a range of q values. Returns a DataFrame where each column is a different q value and each row is a different patient. 

Input: DataFrame of frequencies per patient and range of q 
Output: DataFrame of diversity scores across the range of q per patient 

"""
function calcOverQ(origDF::DataFrame,rangeQ::Array{Float64})
  qDdf = DataFrame(patient=unique(origDF[:,:patient]))

  for q in rangeQ
      tempDF = patientdiversity(origDF,q);
      qDdf[:,Symbol(q)]=tempDF[:,:diversity]
  end
  
  return qDdf
end

"""    findInflection(q::Array, qD::Array)

Computes the inflection point of a sequence numerically.

"""
function findInflection(q::Array, qD::Array)
  approxLogDeriv = q[2:end].*diff(qD)./diff(q);
  return approxLogDeriv; 
end


"""    reportInflection(q::Array, qD::Array)

Computes the inflection point and the magnitude of the slope at that point using findInflection and findInflectionLocal functions.

"""
function reportInflection(q::Array, qD::Array)
  diffVector = findInflection(q, qD)
  inflectPt = argmin(diffVector)
  slopeInflectPt = abs(diffVector[inflectPt])
  qInflecPt = q[inflectPt+1]
  return (qInflecPt, slopeInflectPt)
end


"""    calcIPdiversity(df::DataFrame, qrange::Array)

Computes diversity over a range of q values (qrange) and an input DataFrame of frequencies across types, then numerically determines the inflection point and the magnitude of the slope at that point, and returns a DataFrame with a column for patient name, a column for q at which the inflection point occurs, and a column for the magnitude of the slope at the inflection point. 

Input: DataFrame and array of q values
Output: DataFrame of inflection point q and slope diversity per patient

"""
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

"""    patientdiversity(countframe::DataFrame,q::Real)

Computes a DataFrame of diversity scores across patients from a countframe DataFrame and a single value of q. 

Input: DataFrame and a real number
Output: DataFrame 

"""
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

"""    requirecol(df::DataFrame, colname::AbstractString) 

Raises an assertion error if colname is not a column of DataFrame.

"""
function requirecol(df::DataFrame, colname::AbstractString) 
  @assert colname in names(df) "Expected column named: $colname, names are $(names(df))"
end

rename!(df::DataFrame, s::String, t::Symbol) = rename!(df, Symbol(s)=>t)


"""    countframe(a::CDR3SeqAnalysis, typeDF::DataFrame)

Counts the number of occurances of CDR3 by patients, after checking that patient and CDR3 columns are present in the DataFrame, and returns as a DataFrame of the counts. 

Input: CDR3SeqAnalysis struct and DataFrame of CDR3 recoveries
Output: DataFrame

"""
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

"""    countframefull(a::CDR3SeqAnalysis, df::DataFrame)

Merges df, a DataFrame of counts, with all the physiochemical information associated with the CDR3 sequences. 

Input: CDR3SeqAnalysis struc and DataFrame
Output: DataFrame 

"""
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

"""    summarize(df::CDR3SeqData)

Computes a summary of the CDR3 sequences recovered in the dataset per patient including the total number of recoveries, minimum number on recoveries per sequence, maximum number of recoveries per sequence, average number of recoveries per sequence, median number of recoveries per sequence, and number of unique CDR3 recoveries.

Input: CDR3SeqData struct
Output: DataFrame

"""
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

"""    filteredDataset(::Type{CDR3SeqData}, df::DataFrame, receptorname::String)

Subsets the dataset based on a receptor type of interest

Input: CDR3SeqData type, a DataFrame of counts, and a string that matches the types of the dataset
Output: DataFrame

Example: df_filtered = filteredDataset(CDR3SeqData, df, "TRA")

"""
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

"""    filteredDataset(::Type{CDR3SeqData}, df::DataFrame, receptorname::Vector{String})

Subsets the dataset based on receptor types of interest, allowing users to subset based on combintations of types. 

Input: CDR3SeqData type, a DataFrame of counts, and a string that matches the types of the dataset
Output: DataFrame

Example: df_filtered = filteredDataset(CDR3SeqData, df, ["TRA","TRB"])

"""
function filteredDataset(::Type{CDR3SeqData}, df::DataFrame, receptorname::Vector{String})
  ## subset the type group only 
  recepframe = @from i in df begin
    @where i.Receptor in receptorname
    @select i
    @collect DataFrame
  end;

  println("Number of recoveries: ", size(recepframe))
  println("Cohort size: ", length(unique(recepframe[:,:Filename])))
  
  ds = CDR3SeqData(recepframe)
  
  return ds
end

"""    DiversityScores{T}

Stores a dataframe of diversity scores along with the range of qs that it was computed over.
"""
struct DiversityScores{T} <: Dataset
  qs::T
  data::DataFrame
end

DataFrame(ds::DiversityScores) = ds.data


"""    DiversityScores(countDF::DataFrame, rangeQ::AbstractVector, samplecol=:patient)

Computes point estimates of diversity from a countframe DataFrame, range of q values, and optional sample column identifier. 

Point estimates of diversity of interest: 
* low q diversity (q=0.01)
* high q diversity (q=100)
* delta qD diversity = low q diversity - high q diversity
* inflection point q
* magnitude of slope at the inflection point
* Shannon diversity (log qD at q=1)
* Simpson diversity (1/qD at q=2)

Input: DataFrame of counts, vector of q values, optional sample column symbol (defaults to :patient)
Output: DiversityScores struct

"""
function DiversityScores(countDF::DataFrame, rangeQ::AbstractVector, samplecol=:patient)
  lowDiv = patientdiversity(countDF,rangeQ[1]);
  ids = lowDiv[:, samplecol]; 
  highDiv = patientdiversity(countDF,rangeQ[end]);
  dfIPslope = calcIPdiversity(countDF,rangeQ);
  ShanDiv = patientdiversity(countDF,1.0); 
  SimpDiv = patientdiversity(countDF,2.0); 
  # compile all diversity scores into a single dataframe
  df = DataFrame(patient=ids,
     lowQ=lowDiv[:,:diversity], 
     highQ=highDiv[:,:diversity], 
     deltaqD=(lowDiv[:,:diversity].-highDiv[:,:diversity]), 
     IPq=dfIPslope[:,:InflPt], 
     IPslope=dfIPslope[:,:InflPtS],
     Shan=log.(ShanDiv[:,:diversity]), 
     Simp= 1.0 ./(SimpDiv[:,:diversity]));
  return DiversityScores(rangeQ, df)
end

# function foo(c::SCRNASeqAnalysis, df::DataFrame)
#   #...SCRNASeq code here
# end

end
