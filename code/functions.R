.distmat <- function(x1,y1,x2,y2){
 xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)
 yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
 t(sqrt(xd + yd)) 
}

.update_A <- function(TT, W, W_star, U, V, S, A_curr, P, zeroA_vec, nzA_vec, A_lb_vec, A_ub_vec, L, J){
 Avec <- vec(A_curr)
 C <- 0
 c <- 0
 for (t in 1:(TT-1)){
  Ct_inv <-  solve(S) %x%  (U[[t]] %*% t(U[[t]]))
  C <- C + Ct_inv
  c <- c + Ct_inv %*% vec( solve((U[[t]] %*% t(U[[t]]))) %*% U[[t]] %*%
                            (W_star[[t]] - (W[[t]] + t(V[[t]]) %*% P)))
 }
 
 Vv <- solve(C)
 m <- Vv %*% c
 
 m_cond <- m[-zeroA_vec] -  Vv[-zeroA_vec, zeroA_vec] %*% solve(Vv[zeroA_vec, zeroA_vec]) %*% m[zeroA_vec]
 Vv_cond <- Vv[-zeroA_vec, -zeroA_vec] - Vv[-zeroA_vec, zeroA_vec] %*% solve(Vv[zeroA_vec, zeroA_vec]) %*% Vv[zeroA_vec, -zeroA_vec]
 
 Aprop <- t(rmvnormRcpp(1, t(m_cond), round(Vv_cond, 8)))
 ids <- which(Aprop >= A_lb_vec & Aprop <= A_ub_vec)
 Avec[nzA_vec[ids]] <- Aprop[ids]
 A <- matrix(Avec, nrow = L, ncol = J, byrow = F)
 return(A)
}

.update_P <- function(TT, W, W_star, U, V, S, P_curr, A, zeroP_vec, nzP_vec, P_lb_vec, P_ub_vec, R, J){
 Pvec <- vec(P_curr)
 C <- 0
 c <- 0
 for (t in 1:(TT-1)){
  if(t  >= 1){
   Ct_inv <-  solve(S) %x%  (V[[t]] %*% t(V[[t]]))
   C <- C + Ct_inv
   c <- c + Ct_inv %*% vec( solve((V[[t]] %*% t(V[[t]]))) %*% V[[t]] %*%
                             (W_star[[t]] - (W[[t]] + t(U[[t]]) %*% A)))
  }
 }
 
 Vv <- solve(C)
 m <- Vv %*% c
 
 m_cond <- m[-zeroP_vec] -  Vv[-zeroP_vec, zeroP_vec] %*% solve(Vv[zeroP_vec, zeroP_vec]) %*% m[zeroP_vec]
 Vv_cond <- Vv[-zeroP_vec, -zeroP_vec] - Vv[-zeroP_vec, zeroP_vec] %*% solve(Vv[zeroP_vec, zeroP_vec]) %*% Vv[zeroP_vec, -zeroP_vec]
 
 # Pprop <-  rmvn(1, t(m_cond), round(as.matrix(Vv_cond), 8))
 Pprop <- t(rmvnormRcpp(1, t(m_cond), round(Vv_cond, 8)))
 ids <- which(Pprop >= P_lb_vec & Pprop <= P_ub_vec)
 Pvec[nzP_vec[ids]] <- Pprop[ids]
 P <- matrix(Pvec, nrow = R, ncol = J, byrow = F)
 return(P)
}


.update_delta <- function(pA, J, TT, X, Z, delta_curr, beta, Y, W, s_eta_vec_ls, deg_id_ls,
                          a0_prior = matrix(0, nrow = 1, ncol = 10), a0_sd = 0.5){
 delta <- delta_curr
 for(i in 1:pA){
  for(j in 1:J){
   delta_prop <- delta
   delta_prop[i,j] <- rnorm(1, delta[i,j], 0.2)
   lhr <- 0
   for(t in 1:TT){
    XA_prop <- pnorm(X[[t]] %*% delta_prop)
    XA <- pnorm(X[[t]] %*% delta)
    ZB <- pnorm(Z[[t]] %*% beta)
    deg_ids <- deg_id_ls[[t]]
    lh1 <- sum(log((1- XA_prop)*(1 -ZB)* dlnorm(Y[[t]], W[[t]], s_eta_vec_ls[[t]]))[which(deg_ids == 0)])-
     sum(log((1- XA)*(1- ZB)* dlnorm(Y[[t]], W[[t]], s_eta_vec_ls[[t]]))[which(deg_ids == 0)])
    lh2 <- sum(log((1-XA_prop)*ZB)[which(deg_ids==2)])-
     sum(log((1-XA)*ZB)[which(deg_ids==2)])
    lh3 <- sum(log(XA_prop[which(deg_ids == 1)])) -
     sum(log(XA[which(deg_ids == 1)])) 
    
    lhr <- lhr + lh1 + lh2 + lh3
   }
   if(i == 1){
    lhr <- lhr + dnorm(delta_prop[i,j], a0_prior[i,j],0.5, log = T)-
     dnorm(delta[i,j], a0_prior[i,j], a0_sd, log = T)
   } else{
    lhr <- lhr + dnorm(delta_prop[i,j],a0_sd, log = T)-
     dnorm(delta[i,j], 0, 5, log = T)
   }
   
   if(lhr > log(runif(1))){
    delta <- delta_prop
   }
  }
 }
 return(delta)
}

.update_beta <- function(pB, J, TT, X, Z, delta, beta_curr, Y, W, s_eta_vec_ls, deg_id_ls,
                         b0_prior =matrix(0, nrow = 1, ncol = 10), b0_sd = 1 ){
 beta <- beta_curr
 for(i in 1:pB){
  for(j in 1:J){
   beta_prop <- beta
   beta_prop[i,j] <- rnorm(1, beta[i,j], 0.2)
   lhr <- 0
   for(t in 1:TT){
    XA <- pnorm(X[[t]] %*% delta)
    ZB_prop <- pnorm(Z[[t]] %*% beta_prop)
    ZB <- pnorm(Z[[t]] %*% beta)
    deg_ids <- deg_id_ls[[t]]
    lh1 <- sum(log((1- XA)*(1- ZB_prop)*dlnorm(Y[[t]], W[[t]], s_eta_vec_ls[[t]]))[which(deg_ids == 0)])-
     sum(log((1- XA)*(1- ZB)*dlnorm(Y[[t]], W[[t]], s_eta_vec_ls[[t]]))[which(deg_ids == 0)])
    lh2 <- sum(log((1-XA)*ZB_prop)[which(deg_ids==2)])-
     sum(log((1-XA)*ZB)[which(deg_ids==2)])
    lhr <- lhr + lh1 + lh2
   }
   if( i == 1){
    lhr <- lhr + dnorm(beta_prop[i,j],b0_prior[i,j],b0_sd, log = T)-
     dnorm(beta[i,j], b0_prior[i,j], b0_sd, log = T)
   }else{
    lhr <- lhr + dnorm(beta_prop[i,j],0,5, log = T)-
     dnorm(beta[i,j], 0, 5, log = T)
   }
   
   if(lhr > log(runif(1))){
    beta <- beta_prop
   }
  }
 }
 return(beta)
}

.update_deg <- function(t, ng, J, deg_id_ls, z_ids_ls, X, Z, delta, beta){
 deg_id_t <-  deg_id_ls[[t]]
 z_ids <- z_ids_ls[[t]]
 XA <- pnorm(X[[t]] %*% delta)
 ZB <- pnorm(Z[[t]] %*% beta)
 prob_deg_cond <- rbernoulli(ng*J,((XA)/(XA + (1-XA)*ZB)))
 for(i in z_ids){
  if(rbernoulli(1, prob_deg_cond[i])){
   deg_id_t[i] <- 1
  } else{
   deg_id_t[i] <- 2
  }
 }
 return(deg_id_t)
}


.get_edge_prop <- function(blob.poly.df){
 ng <- nrow(blob.poly.df)
 nb_ls <- Touching_List <- rgeos::gTouches(blob.poly.df, byid = TRUE, returnDense=FALSE)
 perimeters <- sp::SpatialLinesLengths(as(blob.poly.df, "SpatialLines"))
 all.length.list <- lapply(1:length(Touching_List), function(from) {
  lines <- rgeos::gIntersection(blob.poly.df[from,], blob.poly.df[Touching_List[[from]],], byid = TRUE)
  if (class(lines) == "SpatialCollections") {
   list.Lines <- lapply(1:length(Touching_List[[from]]), function(to) {
    line.single <- rgeos::gIntersection(blob.poly.df[from,], blob.poly.df[Touching_List[[from]][to],])
    if (class(line.single) == "SpatialPoints") {
     # Double the point to create a line
     L1 <- rbind(line.single@coords, line.single@coords)
     rownames(L1) <- letters[1:2]
     Sl1 <- Line(L1)
     Lines.single <- Lines(list(Sl1), ID = as.character(to))
    } else if (class(line.single) == "SpatialLines") {
     Lines.single <- line.single@lines[[1]]
     Lines.single@ID <- as.character(to)
    }
    Lines.single
   })
   lines <- SpatialLines(list.Lines)
  }
  l_lines <- sp::SpatialLinesLengths(lines)
  res <- data.frame(origin = from,
                    perimeter = perimeters[from],
                    touching = Touching_List[[from]],
                    t.length = l_lines,
                    t.pc = l_lines/perimeters[from])
  res
 })
 perims<- do.call("rbind", all.length.list) %>%
  filter(t.pc > 0)%>%
  arrange(origin, touching)  %>%
  # dplyr::select(-c(perimeter, t.length))
  dplyr::select(-c(perimeter))
 
 
 perims2 <- left_join(perims %>%
                       rename("blobID" = "origin"), blob.poly.df@data  %>%
                       dplyr::select(-geeID), by = "blobID") %>%
  left_join(.,  blob.poly.df@data %>%
             dplyr::select(-geeID) %>%
             rename("touching" = "blobID", 
                    "lc_touching" = "lc_blob"), by = "touching")
 
 ret <- perims2 %>%
  filter(lc_blob!= lc_touching)%>%
  group_by(blobID)%>%
  summarise(edge_pc = sum(t.pc), edge_length = sum(t.length)) %>%
  left_join(data.frame(blobID= 1:ng),. , by = "blobID") %>%
  mutate(edge_pc = if_else(is.na(edge_pc), 0, edge_pc),
         edge_length = if_else(is.na(edge_length), 0, edge_length),
         edge_length = edge_length/1000)
 return(ret)
}

get_H_mat2 <- function(blob.poly.df, crs.utm =  CRS("+proj=utm +zone=17 +datum=WGS84 +units=km"),
                       fd_prob, df_prob, phi = 1, l= 0.1,gauss_nb = F, n = 500000){
 ng <- nrow(blob.poly.df)
 blob.poly.df <- spTransform(blob.poly.df[match(sort(blob.poly.df$blobID),blob.poly.df$blobID,),],
                             crs.utm) 
 full_id <- data.frame(origin= sort(rep(1:ng, ng)), touching =rep(1:ng, ng))
 lc.list <- lapply(1:ng, function(x){
  res <- blob.poly.df@data %>%
   mutate(startID = x,
          startLC = blob.poly.df@data$lc_blob[x])%>%
   rename("endID" = "blobID", "endLC" = "lc_blob" )%>%
   dplyr::select(startID, endID, startLC, endLC)
  res})

 lc.df <- do.call("rbind", lc.list) %>%
  mutate(lc_prob = case_when(endLC == "developed" ~ fd_prob,
                             endLC == "forest" ~ df_prob))
 # lc:lc transitions
 H_lc <- matrix(lc.df$lc_prob, nrow = ng, ncol = ng, byrow = T)
 
 
 nb_ls <- Touching_List <- rgeos::gTouches(blob.poly.df, byid = TRUE, returnDense=FALSE)
 perimeters <- sp::SpatialLinesLengths(as(blob.poly.df, "SpatialLines"))
 all.length.list <- lapply(1:length(Touching_List), function(from) {
  lines <- rgeos::gIntersection(blob.poly.df[from,], blob.poly.df[Touching_List[[from]],], byid = TRUE)
  if (class(lines) == "SpatialCollections") {
   list.Lines <- lapply(1:length(Touching_List[[from]]), function(to) {
    line.single <- rgeos::gIntersection(blob.poly.df[from,], blob.poly.df[Touching_List[[from]][to],])
    if (class(line.single) == "SpatialPoints") {
     # Double the point to create a line
     L1 <- rbind(line.single@coords, line.single@coords)
     rownames(L1) <- letters[1:2]
     Sl1 <- Line(L1)
     Lines.single <- Lines(list(Sl1), ID = as.character(to))
    } else if (class(line.single) == "SpatialLines") {
     Lines.single <- line.single@lines[[1]]
     Lines.single@ID <- as.character(to)
    }
    Lines.single
   })
   lines <- SpatialLines(list.Lines)
  }
  l_lines <- sp::SpatialLinesLengths(lines)
  res <- data.frame(origin = from,
                    perimeter = perimeters[from],
                    touching = Touching_List[[from]],
                    t.length = l_lines,
                    t.pc = l_lines/perimeters[from])
  res
 })
 all.length.df <- do.call("rbind", all.length.list) %>%
  filter(t.pc > 0)
 
 borders <- all.length.df %>%
  group_by(origin) %>%
  mutate(total = sum(t.pc)) %>%
  ungroup() %>%
  distinct(origin, total, perimeter) %>%
  filter(total  < 1) %>%
  mutate(t.pc = 1-total,
         t.length = t.pc * perimeter, 
         touching = origin) %>%
  dplyr::select(c(origin, perimeter, touching, t.length, t.pc))
 
 perims <- all.length.df %>% 
  arrange(origin, touching)  %>%
  dplyr::select(-c(perimeter, t.length))
 
 perims2 <- left_join(perims %>%
                       rename("blobID" = "origin"), blob.poly.df@data  %>%
                       dplyr::select(-geeID), by = "blobID") %>%
  left_join(.,  blob.poly.df@data %>%
             dplyr::select(-geeID) %>%
             rename("touching" = "blobID", 
                    "lc_touching" = "lc_blob"), by = "touching")
 
 

 area.list <- lapply(1:ng, function(x){
  a <- rgeos::gArea(blob.poly.df[x,])
  res <- data.frame(origin = x,
                    area = a)
  res})
 area.df <- do.call("rbind", area.list) 
 
 
 start_pts <- spsample(blob.poly.df,n=n,"random")
 o_start <- over(start_pts, blob.poly.df)
 temp_pts <- start_pts@coords + rmvnorm(n, rep(0,2), phi*diag(2))
 end_pts <- SpatialPoints(coords = coordinates(temp_pts),
                          proj4string = crs.utm)
 o_end <- over(end_pts, blob.poly.df)
 na_inds <- which(is.na(o_end$blobID))
 count <- 1
 df_ls <- list()
 df_ls[[count]] <-  data.frame(startID = o_start$blobID,
                               endID = o_end$blobID) %>%
  na.omit()
 nn <- nrow(df_ls[[count]])
 while(nn < n){
  start_pts <- spsample(blob.poly.df,n=length(na_inds)+1000,"random")
  o_start <- over(start_pts, blob.poly.df)
  temp_pts <- start_pts@coords + rmvnorm(length(na_inds)+1000, rep(0,2), phi*diag(2))
  end_pts <- SpatialPoints(coords = coordinates(temp_pts),
                           proj4string = crs.utm)
  o_end <- over(end_pts, blob.poly.df)
  na_inds <- which(is.na(o_end$blobID))
  count <- count + 1
  df_ls[[count]] <-  data.frame(startID = o_start$blobID,
                                endID = o_end$blobID) %>%
   na.omit()
  nn <- nn + nrow(df_ls[[count]])
  
 }  
 
 df <- do.call(rbind, df_ls)
 
 
 #rows sum to total originiating from startID
 # row = startID, col = endID
 H_gauss <- as.matrix(table(df$startID, df$endID))
 H_gauss <- H_gauss/rowSums(H_gauss)
 pi<- diag(H_gauss)
 diag(H_gauss) <- 0
 H_local <- H_gauss
 
 # if wanting to restrict gaussian dispersal to neighbors
 if(gauss_nb){
  for(i in 1:ng){
   H_local[i, -c(nb_ls[[i]])] <- 0
  }
 }

 H_local <- H_local*H_lc
 H_local <- H_local/rowSums(H_local)
 
 
 H_long <- matrix(area.df$area, nrow = ng, ncol = ng,
                  byrow = T)
 for (i in 1:ng){
  H_long[i,i] <- 0
  # H_long[i, nb_ls[[i]]] <- 0
 }
 H_long <- H_long * H_lc
 H_long <- H_long/rowSums(H_long)
 
 H <- (diag(pi) +
        matrix((1-pi)*(1-l), nrow = ng, ncol = ng)*H_local+
        matrix((1-pi)*(l), nrow = ng, ncol = ng)*H_long)
 H <- matrix(H,nrow = ng, ncol = ng)
 return(H)

}


