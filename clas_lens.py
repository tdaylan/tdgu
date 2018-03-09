from __init__ import *
  

def test_clas_lens(strgcnfgextnexec=None):
   
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
    indxprev = array(indxprev)
    indxiter = setdiff1d(indxtotl, indxprev)
    indxiter = choice(indxiter, size=numbiter, replace=False)
    
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
                                  listnamecnfgextn, \
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


def retr_conv(inpt, size_in, size_out, name="conv"):
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


def retr_fulc(inpt, size_in, size_out, name="fulc"):
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
        tens = retr_conv(tens, 1, 32, "conv0000")
        tf.summary.image('conv0000', tens[:, :, :, 0:1], 3)
        tens = retr_conv(tens, 32, 64, "conv0001");
        tf.summary.image('conv0001', tens[:, :, :, 0:1], 3)
    else:
        tens = retr_conv(tens, 1, 64, "conv")
        tf.summary.image('conv', tens, 3)
        tens = tf.nn.max_pool(tens, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding="SAME")
        tf.summary.image('convpool', tens, 3)

    tens = tf.reshape(tens, [sizebtch, 25 * 25 * 64])

    if boolfulcdoub:
        tens = retr_fulc(tens, 25 * 25 * 64, sizeembd, "fulc0000")
        inptembd = tens
        logits = retr_fulc(tens, sizeembd, 2, "fulc0001")
    else:
        inptembd = tens
        logits = retr_fulc(tens, 25*25*64, 2, "fulc")

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


def test_line():
    
    import tensorflow as tf
    sizefeat = 100
    truewght = 3.
    truebias = 2.
    pathclaslens = os.environ["TDGU_DATA_PATH"] + '/clas_lens/data/'
    feat = np.linspace(-1, 1, sizefeat)
    truelabl = truewght * feat + truebias + np.random.randn(sizefeat)
    tensfeat = tf.placeholder("float")
    tenstruelabl = tf.placeholder("float")
    tenswght = tf.Variable(np.random.randn(), name="wght")
    tensbias = tf.Variable(np.random.randn(), name="bias")
    tensmodllabl = tf.add(tf.multiply(tensfeat, tenswght), tensbias)
    tf.summary.scalar("wght", tenswght)
    tf.summary.histogram("wght", tenswght)
    tf.summary.scalar("bias", tensbias)
    tenscost = tf.reduce_sum(tf.square(tenstruelabl - tensmodllabl))
    tf.summary.scalar("tenscost", tenscost)
    operloss = tf.train.GradientDescentOptimizer(0.0001).minimize(tenscost)
    smry = tf.summary.merge_all()
    objtwrit = tf.summary.FileWriter(pathclaslens + 'line')
    with tf.Session() as sess:
        tf.global_variables_initializer().run()
        for i in range(1000):
            sess.run(operloss, feed_dict={tensfeat: feat, tenstruelabl: truelabl})
            meta = sess.run(smry, feed_dict={tensfeat: feat, tenstruelabl: truelabl})
            objtwrit.add_summary(meta, i)


def test_conv():

    import os.path
    import shutil
    import tensorflow as tf
    pathclaslens = os.environ["TDGU_DATA_PATH"] + '/clas_lens/data/'
    LABELS = os.path.join(os.getcwd(), "labels_1024.tsv")
    SPRITES = os.path.join(os.getcwd(), "sprite_1024.png")
    mnist = tf.contrib.learn.datasets.mnist.read_data_sets(train_dir=pathclaslens + "data", one_hot=True)
   
    ratelern = 0.001
    tf.reset_default_graph()
    sess = tf.Session()
    
    x = tf.placeholder(tf.float32, shape=[None, 784], name="x")
    y = tf.placeholder(tf.float32, shape=[None, 10], name="truelabl")
    feat = tf.reshape(x, [-1, 28, 28, 1])
    
    tf.summary.image('input', feat, 3)
    
    tens = retr_conv(feat, 1, 32, "conv0000")
    tens = retr_conv(tens, 32, 64, "conv0001")
    tens = tf.reshape(tens, [-1, 7 * 7 * 64])
    fc1 = retr_fulc(tens, 7 * 7 * 64, 1024, "fc1")
    relu = tf.nn.relu(fc1)
    embedding_input = relu
    tf.summary.histogram("fc1/relu", relu)
    sizeembd = 1024
    logits = retr_fulc(fc1, 1024, 10, "fc2")
    
    with tf.name_scope("cent"):
      cent = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=logits, labels=y), name="cent")
      tf.summary.scalar("cent", cent)
    
    with tf.name_scope("train"):
      opertran = tf.train.AdamOptimizer(ratelern).minimize(cent)
    
    with tf.name_scope("accu"):
      accu = tf.reduce_mean(tf.cast(tf.equal(tf.argmax(logits, 1), tf.argmax(y, 1)), tf.float32))
      tf.summary.scalar("accu", accu)
    
    summ = tf.summary.merge_all()
    
    embedding = tf.Variable(tf.zeros([1024, sizeembd]), name="test_embedding")
    assignment = embedding.assign(embedding_input)
    saver = tf.train.Saver()
    
    sess.run(tf.global_variables_initializer())
    writer = tf.summary.FileWriter(pathclaslens + 'conv')
    writer.add_graph(sess.graph)
    
    cnfgtdbr = tf.contrib.tensorboard.plugins.projector.ProjectorConfig()
    cnfgembd = cnfgtdbr.embeddings.add()
    cnfgembd.tensor_name = embedding.name
    cnfgembd.sprite.image_path = SPRITES
    cnfgembd.metadata_path = LABELS
    cnfgembd.sprite.single_image_dim.extend([28, 28])
    tf.contrib.tensorboard.plugins.projector.visualize_embeddings(writer, cnfgtdbr)
    for i in range(200):
        batch = mnist.train.next_batch(100)
        if i % 5 == 0:
            tensaccu, s = sess.run([accu, summ], feed_dict={x: batch[0], y: batch[1]})
            writer.add_summary(s, i)
        #if i % 500 == 0:
        #  sess.run(assignment, feed_dict={x: mnist.test.images[:1024], y: mnist.test.labels[:1024]})
        #  saver.save(sess, os.path.join(pathclaslens, "model.ckpt"), i)
        sess.run(opertran, feed_dict={x: batch[0], y: batch[1]})


def test_mlor():
    
    import tensorflow as tf
    
    tf.reset_default_graph()
    
    sizebthc = 100
    ratelern = 0.5
    numbepoc = 5
    pathclaslens = os.environ["TDGU_DATA_PATH"] + '/clas_lens/data/'
    
    from tensorflow.examples.tutorials.mnist import input_data
    mnist = input_data.read_data_sets('MNIST_data', one_hot=True)
    
    with tf.name_scope('inpt'):
        tensfeat = tf.placeholder(tf.float32, shape=[None, 784], name="feat") 
        tenstruelabl = tf.placeholder(tf.float32, shape=[None, 10], name="truelabl")
    
    with tf.name_scope("wght"):
        variwght = tf.Variable(tf.zeros([784, 10]))
    
    with tf.name_scope("bias"):
        varibias = tf.Variable(tf.zeros([10]))
    
    with tf.name_scope("softmax"):
        tensmodllabl = tf.nn.softmax(tf.matmul(tensfeat, variwght) + varibias)
    
    with tf.name_scope('cent'):
        tenscent = tf.reduce_mean(-tf.reduce_sum(tenstruelabl * tf.log(tensmodllabl), reduction_indices=[1]))
    
    with tf.name_scope('tran'):
        opertran = tf.train.GradientDescentOptimizer(ratelern).minimize(tenscent)
    
    with tf.name_scope('accu'):
        tensaccu = tf.reduce_mean(tf.cast(tf.equal(tf.argmax(tensmodllabl, 1), tf.argmax(tenstruelabl, 1)), tf.float32))
        
    tf.summary.scalar("cent", tenscent)
    tf.summary.scalar("accu", tensaccu)
    
    opersmry = tf.summary.merge_all()
    
    numbbtch = int(mnist.train.num_examples / sizebthc)
    
    with tf.Session() as sess:
        sess.run(tf.global_variables_initializer())
        objtwrit = tf.summary.FileWriter(pathclaslens + 'mlor')
        objtwrit.add_graph(sess.graph)
        for k in range(numbepoc):
            for i in range(numbbtch):
                featbtch, truelablbtch = mnist.train.next_batch(sizebthc)
                temp, summary = sess.run([opertran, opersmry], feed_dict={tensfeat: featbtch, tenstruelabl: truelablbtch})
                objtwrit.add_summary(summary, k * numbbtch + i)


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
