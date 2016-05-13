#See instruction  https://www.youtube.com/watch?v=lT4Kosc_ers

#install.packages("twitteR")
#install.packages("RCurl")
#install.packages("tm")
#install.packages("wordcloud")

require(twitteR)     # Twitter API (Application Program Interface)
require(RCurl)       # Accessing internet from R
require(tm)          # text mining
require(wordcloud)   # wordcloud

# create API creditial through twitter at http://apps.twitter.com
# Mobile phone number is needed: Settings-->Your Twitter Data--> Phone
# Mind_of_Cities"
consumer_key <- "RXXXXXXXYUa5Fpzl"
consumer_secret <- "68dXXXXXXXmhQsOG9"
access_token <-  "34157XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXe0AwoT3QoiQ40P"
access_secret <- "Tc7ASXXXXXXXXXXXXXXXXXXXXXXXXtHVULPdD5"

#establish connection
setup_twitter_oauth(consumer_key, consumer_secret, access_token, access_secret)

# this function remove non-ACCII characters from a string
Remove_nonASCII <- function (x) {
# convert string to vector of words
x <- unlist(strsplit(x, split=" "))
# find indices of words with non-ASCII characters
nonAscIDX <- grep("x", iconv(x, "latin1", "ASCII", sub="x"))
# subset original vector of words to exclude words with non-ACCII characters
if( length(nonAscIDX) >0)
  ascVec <- x[ - nonAscIDX]  else 
 ascVec <- x;
# convert vector back to string
x <- paste(ascVec, collapse = " ")
return(x)
}

# this function cleans up a body of text
clean_text <- function (docs) {
#remove non-English characters (ASCII)
docs = unlist( lapply( docs,Remove_nonASCII))
#clean up text according to
# http://www.sthda.com/english/wiki/text-mining-and-word-cloud-fundamentals-in-r-5-simple-steps-you-should-know
docs = Corpus(VectorSource(docs))
# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)
# Text stemming
# docs <- tm_map(docs, stemDocument)
return(docs)
}

tweets = searchTwitter("Clinton",lang="en",n=500, resultType="recent")
# extract text from tweets
texts = sapply(tweets, function(x) x$getText())
docs = clean_text(texts)
# Remove your  query word
docs <- tm_map(docs, removeWords, c("hillary", "clinton")) 
# default wordcloud
wordcloud(docs)
# fine-tuned
wordcloud(docs, random.order=F, max.words =40,col=rainbow(40))


# what's on the minds of poeple in Sioux Falls, within 5 mile ? 
# geocode lookup:  http://mygeoposition.com/
tweets = searchTwitter(" ",lang="en",n=500, geocode='43.544596,-96.731103,5mi', resultType="recent")
docs = clean_text(tweets)
wordcloud(docs)

# other ways to search twitter   https://dev.twitter.com/rest/public/search
tweets = searchTwitter("#beer", n=100)  #BoG16
## Search between two dates
tweets =searchTwitter('breast+cancer', since='2016-05-09', until='2016-05-13')

## using resultType
tweets =searchTwitter('world cup+brazil', resultType="popular", n=15)
tweets =searchTwitter('from:realDonaldTrump', resultType="recent", n=500)



