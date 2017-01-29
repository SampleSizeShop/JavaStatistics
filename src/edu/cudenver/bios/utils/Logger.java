/*
 * Java Statistics.  A java library providing power/sample size estimation for
 * the general linear model.
 *
 * Copyright (C) 2017 Regents of the University of Colorado.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
package edu.cudenver.bios.utils;

//import java.util.function.Supplier; // in Java 8
import org.apache.log4j.Level;

/**
 * This class wraps the Log4J Logger class to provide methods
 * that take "message suppliers" rather than messages. A message
 * supplier is a Supplier that supplies messages on demand.
 * Thus the evaluation of messages that are expensive to evaluate
 * can be deferred until after determining that the logger is
 * enabled for the level.
 *
 * <p>
 * As a convenience, this class also provides "pass-through"
 * methods that simply delegate to the wrapped logger.
 */
public class Logger {
    /**
     * The wrapped logger: the logger to which logging method
     * invocations are delegated.
     */
    private final org.apache.log4j.Logger logger;

    /**
     * Retrieve an instance of this class wrapping a logger
     * named according to the name of a class.
     *
     * @param c The class.
     *
     * @return The instance.
     */
    public static Logger getLogger(Class c) {
        return new Logger(org.apache.log4j.Logger.getLogger(c));
    }

    /**
     * Construct an instance of this class.
     *
     * @param logger The wrapped logger.
     */
    private Logger(org.apache.log4j.Logger logger) {
        this.logger = logger;
    }

    /**
     * Log a message object with the DEBUG level.
     *
     * @param message The message object to log.
     */
    public void debug(Object message) {
        logger.debug(message);
    }

    /**
     * Log a message object with the DEBUG level including the stack trace
     * of the Throwable t passed as parameter.
     *
     * @param message The message object to log.
     * @param t       The exception to log, including its stack trace.
     */
    public void debug(Object message, Throwable t) {
        logger.debug(message, t);
    }

    /**
     * Log a supplied message object with the DEBUG level.
     *
     * @param supplier The supplier of the message object to log.
     */
    public void debug(Supplier<Object> supplier) {
        if (logger.isDebugEnabled()) {
            logger.debug(supplier == null ? null : supplier.get());
        }
    }

    /**
     * Log a supplied message object with the DEBUG level including the stack trace
     * of the Throwable t passed as parameter.
     *
     * @param supplier The supplier of the message object to log.
     * @param t        The exception to log, including its stack trace.
     */
    public void debug(Supplier<Object> supplier, Throwable t) {
        if (logger.isDebugEnabled()) {
            logger.debug(supplier == null ? null : supplier.get(), t);
        }
    }

    /**
     * Log a message object with the INFO level.
     *
     * @param message The message object to log.
     */
    public void info(Object message) {
        logger.info(message);
    }

    /**
     * Log a message object with the INFO level including the stack trace
     * of the Throwable t passed as parameter.
     *
     * @param message The message object to log.
     * @param t       The exception to log, including its stack trace.
     */
    public void info(Object message, Throwable t) {
        logger.info(message, t);
    }

    /**
     * Log a supplied message object with the INFO level.
     *
     * @param supplier The supplier of the message object to log.
     */
    public void info(Supplier<Object> supplier) {
        if (logger.isInfoEnabled()) {
            logger.info(supplier == null ? null : supplier.get());
        }
    }

    /**
     * Log a supplied message object with the INFO level including the stack trace
     * of the Throwable t passed as parameter.
     *
     * @param supplier The supplier of the message object to log.
     * @param t        The exception to log, including its stack trace.
     */
    public void info(Supplier<Object> supplier, Throwable t) {
        if (logger.isInfoEnabled()) {
            logger.info(supplier == null ? null : supplier.get(), t);
        }
    }

    /**
     * Log a message object with the WARN level.
     *
     * @param message The message object to log.
     */
    public void warn(Object message) {
        logger.warn(message);
    }

    /**
     * Log a message object with the WARN level including the stack trace
     * of the Throwable t passed as parameter.
     *
     * @param message The message object to log.
     * @param t       The exception to log, including its stack trace.
     */
    public void warn(Object message, Throwable t) {
        logger.warn(message, t);
    }

    /**
     * Log a supplied message object with the WARN level.
     *
     * @param supplier The supplier of the message object to log.
     */
    public void warn(Supplier<Object> supplier) {
        if (logger.isEnabledFor(Level.WARN)) {
            logger.warn(supplier == null ? null : supplier.get());
        }
    }

    /**
     * Log a supplied message object with the WARN level including the stack trace
     * of the Throwable t passed as parameter.
     *
     * @param supplier The supplier of the message object to log.
     * @param t        The exception to log, including its stack trace.
     */
    public void warn(Supplier<Object> supplier, Throwable t) {
        if (logger.isEnabledFor(Level.WARN)) {
            logger.warn(supplier == null ? null : supplier.get(), t);
        }
    }
}
