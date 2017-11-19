from __init__ import *
  

def pcat_clas_lens(strgcnfgextnexec=None):
   
    pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
    
    dictargs = {}
    dictargs['elemtype'] = ['lens']
    dictargs['exprtype'] = 'hubb'
    dictargs['elemtype'] = ['lens']
    dictargs['maxmnumbelempop0reg0'] = 0
    dictargs['numbelempop0reg0'] = 0
    dictargs['makeplot'] = False
    dictargs['mockonly'] = True
    dictargs['verbtype'] = 0
    
    dictargsvari = {}
    numbiter = 1000
    numbtotl = 1000000
    indxtotl = arange(numbtotl)

    indxprev = []
    listrtag = fnmatch.filter(os.listdir(pathimag), '*pcat_clas_lens_*')
    for rtag in listrtag:
        boolgdatinit = pcat.util.chec_statfile(rtag, 'gdatinit')
        if not boolgdatinit:
            pcat.util.dele_rtag(rtag)
        else:
            indxprev.append(int(rtag.split('pcat_clas_lens_')[1][4:12]))
    print 'indxtotl'
    summgene(indxtotl)
    indxprev = array(indxprev)
    print 'indxprev'
    summgene(indxprev)
    indxiter = setdiff1d(indxtotl, indxprev)
    print 'indxiter'
    summgene(indxiter)
    indxiter = choice(indxiter, size=numbiter, replace=False)
    print 'indxiter'
    summgene(indxiter)
            
    
    for k in indxiter:
        if rand() > 0.5:
            namecnfgextn =  'none%08d' % k
        else:
            namecnfgextn =  'lens%08d' % k
        
        dictargsvari[namecnfgextn] = {}
        if namecnfgextn.startswith('lens'):
            dictargsvari[namecnfgextn]['truefluxsourreg0'] = 1e-22
        dictargsvari[namecnfgextn]['seedtype'] = k
    
    dictglob = pcat.main.initarry( \
                                  dictargsvari, \
                                  dictargs, \
                                  strgcnfgextnexec=strgcnfgextnexec, \
                                 )


def writ_data():
    import tensorflow as tf
    
    pathdata = os.environ["PCAT_DATA_PATH"] + '/data/outp/'
    listrtagnone = fnmatch.filter(os.listdir(pathdata), '20*pcat_clas_lens_none*')
    listrtaglens = fnmatch.filter(os.listdir(pathdata), '20*pcat_clas_lens_lens*')
    listrtag = listrtagnone + listrtaglens
    boollenstemp = []
    cntpdatatemp = []
    for k, rtag in enumerate(listrtag):
        print 'Processing %s...' % rtag
        
        boolgdatinit = pcat.util.chec_statfile(rtag, 'gdatinit')
        if not boolgdatinit:
            continue
        
        pathoutprtag = pcat.util.retr_pathoutprtag(rtag)
        path = pathoutprtag + 'gdatinit'
        gdat = pcat.util.readfile(path) 
        cntpdatatemp.append(gdat.cntpdatareg0[0, :, 0])
        if rtag in listrtaglens:
            boollenstemp.append(True)
        else:
            boollenstemp.append(False)
    
    numbcnfg = len(boollenstemp)
    cntpdata = empty((numbcnfg, cntpdatatemp[0].size))
    boollens = zeros(numbcnfg, dtype=bool)
    indxcnfg = arange(numbcnfg)
    print 'numbcnfg'
    print numbcnfg
    for k in indxcnfg:
        boollens[k] = boollenstemp[k]
        cntpdata[k, :] = cntpdatatemp[k]
        print 'k'
        print k
        print 'cntpdata[k, :]'
        summgene(cntpdata[k, :])
        print 'boollens[k]'
        print boollens[k]
        print

    pathclaslens = os.environ["TDGU_DATA_PATH"] + '/clas_lens/data/'
    path = pathclaslens + 'claslens.h5'
    print 'Writing to %s...' % path
    filearry = h5py.File(path, 'w')
    filearry.create_dataset('cntpdata', data=cntpdata)
    filearry.create_dataset('boollens', data=boollens)
    filearry.close()


def conv_layer(inpt, size_in, size_out, name="conv"):
  import tensorflow as tf
  with tf.name_scope(name):
      wght = tf.Variable(tf.truncated_normal([5, 5, size_in, size_out], stddev=0.1), name="wght" + name)
      bias = tf.Variable(tf.constant(0.1, shape=[size_out]), name="bias" + name)
      conv = tf.nn.conv2d(inpt, wght, strides=[1, 1, 1, 1], padding="SAME")
      acti = tf.nn.relu(conv + bias)
      tf.summary.histogram("wght" + name, wght)
      tf.summary.histogram("bias" + name, bias)
      tf.summary.histogram("acti" + name, acti)
      return tf.nn.max_pool(acti, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding="SAME")


def fc_layer(inpt, size_in, size_out, name="fulc"):
  import tensorflow as tf
  with tf.name_scope(name):
      wght = tf.Variable(tf.truncated_normal([size_in, size_out], stddev=0.1), name="wght" + name)
      bias = tf.Variable(tf.constant(0.1, shape=[size_out]), name="bias" + name)
      acti = tf.nn.relu(tf.matmul(inpt, wght) + bias)
      tf.summary.histogram("wght" + name, wght)
      tf.summary.histogram("bias" + name, bias)
      tf.summary.histogram("acti" + name, acti)
      return acti


def clas_lens_wrap(ratelern, boolconvdoub, boolfulcdoub, strghypr):
    
    print 'Classifer initialized.'

    import tensorflow as tf

    boolexeclens = False

    if boolexeclens:
        pathclaslens = os.environ["TDGU_DATA_PATH"] + '/clas_lens/data/'
        path = pathclaslens + 'claslens.h5'
        print 'Reading %s...' % path
        filearry = h5py.File(path, 'r')
        cntpdata = filearry['cntpdata'][()]
        boollens = filearry['boollens'][()]
        filearry.close()
   
        numbdata = cntpdata.shape[0]
    else:
        from tensorflow.examples.tutorials.mnist import input_data
        mnist = input_data.read_data_sets('MNIST_data', one_hot=True)
        cntpdata = mnist[0]
        boollens = mnist[1]
    
    print 'Found %s images.' % numbdata
    indxdata = arange(numbdata)

    if boolexeclens:
        print 'Randomizing the order of the dataset...'
        indxrndm = choice(indxdata, size=numbdata, replace=False)
        cntpdata = cntpdata[indxrndm, :]
        boollens = boollens[indxrndm]

        # make labels one-hot encoded
        boollenstemp = copy(boollens)
        boollens = zeros((numbdata, 2))
        boollens[where(logical_not(boollenstemp))[0], 0] = 1.
        boollens[where(boollenstemp)[0], 1] = 1.


    numbdatatran = int(0.5 * numbdata) - int(0.5 * numbdata) % 10
    sizebtch = numbdatatran / 10
    
    print 'Using %d of the images as the training dataset.' % numbdatatran 
    print 'Will train in 10 batches with batch size %d' % sizebtch
    numbside = int(sqrt(cntpdata.shape[1]))
    numbdatatest = numbdata - numbdatatran
    indxdatatest = arange(numbdatatest)
    filelabl = open(pathclaslens + 'labl.tsv', 'w')
    for k in indxdatatest:
        if boollens[k, 0] == 1.:
            filelabl.write('0\t')
        if boollens[k, 1] == 1.:
            filelabl.write('1\t')
    filelabl.close()
    
    numbimagsidesprt = int(ceil(sqrt(numbdatatest)))
    numbsidesprt = numbside * numbimagsidesprt
    cntpdatasprt = zeros((numbsidesprt, numbsidesprt))
    for k in indxdatatest:
        indxxaxi = k % numbimagsidesprt
        indxyaxi = k // numbimagsidesprt
        cntpdatasprt[indxyaxi*numbside:(indxyaxi+1)*numbside, indxxaxi*numbside:(indxxaxi+1)*numbside] = cntpdata[k, :].reshape((numbside, numbside))
        
    cntpdatasprt /= amax(cntpdatasprt)
    cntpdatasprt *= 255.
    sp.misc.imsave(pathclaslens + 'sprt.png', cntpdatasprt)
     
    sizeembd = 1024

    tf.reset_default_graph()
    sess = tf.Session()

    tenscntpdata = tf.placeholder(tf.float32, shape=[sizebtch, 10000], name="cntpdata")
    tenstruelabl = tf.placeholder(tf.float32, shape=[sizebtch, 2], name="labl")
    
    tens = tf.reshape(tenscntpdata, [sizebtch, 100, 100, 1])
    tf.summary.image('inpt', tens, 3)

    if boolconvdoub:
        tens = conv_layer(tens, 1, 32, "conv0000")
        tf.summary.image('conv0000', tens[:, :, :, 0:1], 3)
        tens = conv_layer(tens, 32, 64, "conv0001");
        tf.summary.image('conv0001', tens[:, :, :, 0:1], 3)
    else:
        tens = conv_layer(tens, 1, 64, "conv")
        tf.summary.image('conv', tens, 3)
        tens = tf.nn.max_pool(tens, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding="SAME")
        tf.summary.image('convpool', tens, 3)

    tens = tf.reshape(tens, [sizebtch, 25 * 25 * 64])

    if boolfulcdoub:
        tens = fc_layer(tens, 25 * 25 * 64, sizeembd, "fulc0000")
        inptembd = tens
        logits = fc_layer(tens, sizeembd, 2, "fulc0001")
    else:
        inptembd = tens
        logits = fc_layer(tens, 25*25*64, 2, "fulc")

    tens = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=logits, labels=tenstruelabl), name="cent")
    tf.summary.scalar("tenscent", tens)

    funcopti = tf.train.AdamOptimizer(ratelern).minimize(tens)

    correct_prediction = tf.equal(tf.argmax(logits, 1), tf.argmax(tenstruelabl, 1))
    accu = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    tf.summary.scalar("tensaccu", accu)

    smry = tf.summary.merge_all()

    embedding = tf.Variable(tf.zeros([sizebtch, sizeembd]), name="embd")
    assignment = embedding.assign(inptembd)
    objtsave = tf.train.Saver()

    sess.run(tf.global_variables_initializer())
    objtwrit = tf.summary.FileWriter(pathclaslens + strghypr)
    objtwrit.add_graph(sess.graph)

    cnfgtsbr = tf.contrib.tensorboard.plugins.projector.ProjectorConfig()
    cnfgembd = cnfgtsbr.embeddings.add()
    cnfgembd.tensor_name = embedding.name
    cnfgembd.sprite.image_path = pathclaslens + 'sprt.png'
    cnfgembd.metadata_path = pathclaslens + 'labl.tsv'
    cnfgembd.sprite.single_image_dim.extend([100, 100])
    tf.contrib.tensorboard.plugins.projector.visualize_embeddings(objtwrit, cnfgtsbr)

    numbiter = numbdatatran / sizebtch
    indxiter = arange(numbiter)
    numbepoc = 5
    indxepoc = arange(numbepoc)
    cntr = 0
    for k in indxepoc:
        for i in indxiter:
            
            indxcnfgbtch = arange(i*sizebtch, (i+1)*sizebtch) + numbdatatest
            batch = [cntpdata[indxcnfgbtch, :], boollens[indxcnfgbtch, :]]
            if cntr % 1 == 0:
                temp, meta = sess.run([accu, smry], feed_dict={tenscntpdata: batch[0], tenstruelabl: batch[1]})
                objtwrit.add_summary(meta, cntr)
            if False and cntr % 1 == 0:
                sess.run(assignment, feed_dict={tenscntpdata: cntpdata[:sizebtch, :], tenstruelabl: boollens[:sizebtch, :]})
                objtsave.save(sess, os.path.join(pathclaslens, "model.ckpt"), cntr)
            print 'cntr: ', cntr
            sess.run(funcopti, feed_dict={tenscntpdata: batch[0], tenstruelabl: batch[1]})
            cntr += 1


def clas_lens():

    for ratelern in [1e-4]:
        for boolfulcdoub in [True]:
            for boolconvdoub in [True]:
                if boolconvdoub:
                    strgconv = 'numbconv0002'
                else:
                    strgconv = 'numbconv0001'
                if boolfulcdoub:
                    strgfulc = 'numbfulc0002'
                else:
                    strgfulc = 'numbfulc0001'
                strghypr = 'ratelern%04g%s%s' % (ratelern, strgconv, strgfulc)
                
                print 'Starting run for %s...' % strghypr
                clas_lens_wrap(ratelern, boolfulcdoub, boolconvdoub, strghypr)


globals().get(sys.argv[1])(*sys.argv[2:])
